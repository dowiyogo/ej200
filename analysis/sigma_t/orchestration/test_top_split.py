"""Meaningful leakage/order/zero-event regression checks; no optical simulation."""
import copy
import unittest
import numpy as np
import top_split as top


def sample():
    rng = np.random.default_rng(4021)
    events, gids, times = [], [], []
    for event in range(160):
        if event == 159: continue  # generated, completely invisible event
        channels = range(16,20) if event % 2 == 0 else range(16,24)
        for channel in channels:
            for i in range(8):
                events.append(event); gids.append(channel)
                times.append(1.0+0.01*i+rng.normal(0,0.03))
    return dict(event_id=np.array(events), global_id=np.array(gids),
                face_type=np.full(len(events),2), time_ns=np.array(times))


class SplitTests(unittest.TestCase):
    def setUp(self):
        self.arr = sample()

    def test_root_order_does_not_affect_channels_or_timestamps(self):
        order = np.random.default_rng(5).permutation(len(self.arr['event_id']))
        shuffled = {k: v[order] for k,v in self.arr.items()}
        self.assertEqual(top.select_channels(self.arr, 0), top.select_channels(shuffled,0))
        for n in (1, 7, 20):
            np.testing.assert_array_equal(top.timestamps(self.arr,[16,17,18,19],0,160,n),
                                          top.timestamps(shuffled,[16,17,18,19],0,160,n))

    def test_eval_poison_cannot_change_training_selection_or_state(self):
        cfg = dict(top.CONFIG, N_BOOTSTRAP=0)
        before = top.learn(self.arr,160,0,cfg)
        poisoned = copy.deepcopy(self.arr)
        mask = poisoned['event_id'] % 2 == 1
        poisoned['global_id'][mask] = 85
        poisoned['time_ns'][mask] += 1000
        after = top.learn(poisoned,160,0,cfg)
        self.assertEqual(before['channels'],after['channels'])
        self.assertEqual(before['winner_N'],after['winner_N'])
        self.assertEqual(before['states'],after['states'])
        prior = copy.deepcopy(before['states'])
        evaluated = top.evaluate(poisoned,160,before)
        self.assertEqual(before['states'], prior)
        self.assertTrue(all(r['n_eff']==0 for r in evaluated['curve']))

    def test_explicit_generated_denominator(self):
        values = top.timestamps(self.arr,[16,17,18,19],1,160,1)
        self.assertEqual(len(values),79)
        self.assertEqual(len(top.generated_ids(160,1)),80)

    def test_eval_cannot_relearn_fit_seeds_or_histogram_axis(self):
        values = top.timestamps(self.arr,[16,17,18,19],0,160,1)
        state = top.fit_state(values, top.CONFIG)
        original = top.fit_engine.gather_seeds
        with top.frozen_fit_construction(state):
            self.assertEqual(top.fit_engine.gather_seeds(values+100),state['seeds'])
            h = top.fit_engine.ROOT.TH1F('test_frozen_axis','',999,-500,500)
            self.assertEqual(h.GetNbinsX(),state['seeds']['n_bins'])
            self.assertEqual(h.GetXaxis().GetXmin(),state['histogram_lo'])
            h.SetDirectory(0)
        self.assertIs(top.fit_engine.gather_seeds, original)


if __name__ == '__main__': unittest.main()

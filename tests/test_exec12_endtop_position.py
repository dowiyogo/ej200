import importlib.util,pathlib,numpy as np
P=pathlib.Path(__file__).resolve().parents[1]/"analysis/exec12_endtop_position.py";S=importlib.util.spec_from_file_location("e12",P);M=importlib.util.module_from_spec(S);S.loader.exec_module(M)
def test_geometry_maps():
 assert np.allclose([M.end_y(i) for i in range(8)],[-26.25,-18.75,-11.25,-3.75,3.75,11.25,18.75,26.25])
 assert M.top_x(0)==-692 and M.top_x(34)==-12 and M.top_x(35)==12 and M.top_x(69)==692
 assert M.top_x(35)-M.top_x(34)==24
def test_pooled_fourth_and_missing():
 assert M.pooled_fourth(np.array([0,0,0,0,1]),np.array([4.,1.,3.,2.,8.]),2)[0]==4
 assert np.isnan(M.pooled_fourth(np.array([0,0,0]),np.array([1.,2.,3.]),1)[0])
def test_weighted_time_centroids():
 t=np.array([[1.,3.,np.nan]]);c=np.array([[4,8,10]]);assert np.isclose(M.weighted_time(t,c)[0],(4+24)/12)
 assert np.isclose(M.centroid(np.array([[1,3]]),np.array([-1.,1.]))[0],.5)
def test_blue_and_singular():
 assert np.allclose(M.blue(np.eye(2)),[.5,.5])
 with np.testing.assert_raises(np.linalg.LinAlgError):M.blue(np.ones((2,2)))
def test_local_pair_gap_and_edge():
 c=np.zeros((2,70),int);t=np.full((2,70),np.nan);c[0,34]=10;c[0,35]=9;c[1,0]=10;c[1,1]=9
 _,_,s,ids=M.local_pair(c,t);assert s[0]=="gap_ambiguous";assert tuple(ids[1])==(16,17)
def test_leave_one_position_concept():
 positions=np.array([-1,0,1]); test=1; train=np.delete(positions,test);assert positions[test] not in train and len(train)==2

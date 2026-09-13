"""Use preserved END functions; parameterize only constants for sensitivity."""
from functools import lru_cache
import json
from pathlib import Path

SOURCE = Path(__file__).resolve().parent/'../upstream/related/420addf/analysis/congruent_sum4_timing.C'


def sensitivity_source():
    source = SOURCE.read_text()
    # A verbatim contiguous slice contains the four complete function bodies.
    begin = source.index('double SprPeakTime()')
    end = source.index('int GroupIndex(',begin)
    return source[begin:end]


@lru_cache(maxsize=1)
def load():
    import ROOT
    ROOT.gROOT.SetBatch(True)
    assert ROOT.gInterpreter.Declare('#include '+json.dumps(str(SOURCE.resolve())))
    assert ROOT.gInterpreter.Declare('namespace exec36param { double kSprRiseNs=.5, kSprFallNs=5., kThresholdPe=4.;\n'+sensitivity_source()+'\n}')
    assert ROOT.gInterpreter.Declare(r'''
    #include <stdexcept>
    namespace exec36 {
    std::vector<double> fit(const std::vector<double>& values, const std::string& name) {
        const auto r=FitCore(values,name);
        return {r.sigmaPs,r.sigmaErrPs,r.chi2Ndf,double(r.n),double(r.usedFit),r.meanNs};
    }
    double original(const std::vector<double>& v) { return LeadingEdgeTime(v); }
    double varied(const std::vector<double>& v,double rise,double fall,double threshold) {
        exec36param::kSprRiseNs=rise;exec36param::kSprFallNs=fall;exec36param::kThresholdPe=threshold;
        return exec36param::LeadingEdgeTime(v);
    }
    int group(int id) {return GroupIndex(id);}
    struct Data {
        std::vector<std::array<std::vector<double>,16>> events;
        long long left=0,right=0;
        Data(const std::string& path,const std::string& cache="") : events(10000) {
            TFile input(path.c_str(),"READ");
            if(input.IsZombie()) throw std::runtime_error("Unreadable input ROOT");
            auto* tree=dynamic_cast<TTree*>(input.Get("end_hits"));
            if(!tree) tree=dynamic_cast<TTree*>(input.Get("sipm_hits"));
            if(!tree) throw std::runtime_error("Missing native/cache hit tree");
            TTreeReader reader(tree);
            TTreeReaderValue<int> event(reader,"event_id"), gid(reader,"global_id"), face(reader,"face_type");
            TTreeReaderValue<double> time(reader,"time_ns");
            std::unique_ptr<TFile> output;
            std::unique_ptr<TTree> cached;
            int ev,id,side;double arrival;
            if(!cache.empty()) {
                output.reset(TFile::Open(cache.c_str(),"CREATE"));
                if(!output || output->IsZombie()) throw std::runtime_error("Cannot create new cache");
                cached.reset(new TTree("end_hits","Exact existing native END hits; no simulation"));
                cached->SetDirectory(nullptr);
                cached->Branch("event_id",&ev);cached->Branch("global_id",&id);
                cached->Branch("face_type",&side);cached->Branch("time_ns",&arrival);
                // Associate baskets with output to avoid keeping the entire tree in RAM.
                cached->SetDirectory(output.get());
            }
            while(reader.Next()) {
                if(*face==2) continue;
                if(*event<0 || *event>=10000 || *gid<0 || *gid>=16 || !std::isfinite(*time)
                   || (*face!=0 && *face!=1) || (*gid<8)!=(*face==0))
                    throw std::runtime_error("END input event/face/ID/time mismatch");
                events[*event][*gid].push_back(*time);
                if(*gid<8) ++left; else ++right;
                if(cached) {ev=*event;id=*gid;side=*face;arrival=*time;cached->Fill();}
            }
            if(cached) {output->cd();cached->Write();cached->SetDirectory(nullptr);output->Close();}
        }
        std::vector<double> times(const std::vector<int>& ids,double rise=.5,double fall=5.,double threshold=4.,bool sensitivity=false) {
            std::vector<double> result;result.reserve(events.size());
            for(const auto& event:events) {
                std::vector<double> arrivals;
                for(int id:ids) {if(id<0 || id>=16) throw std::runtime_error("Bad group ID");
                    arrivals.insert(arrivals.end(),event[id].begin(),event[id].end());}
                const double value=sensitivity?varied(arrivals,rise,fall,threshold):LeadingEdgeTime(arrivals);
                result.push_back(value);
            }
            return result;
        }
    };
    }
    ''')
    return ROOT

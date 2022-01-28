#include "specToChrom.h"
saristc/mzMLTrace
#include <iostream>
saristc/mzMLTrace

saristc/mzMLTrace
namespace asaristc {
saristc/mzMLTrace

saristc/mzMLTrace
chromatogram::chromatogram(int &count, point &peak, double * mzs, double * rts, double * intns) {
saristc/mzMLTrace
  mz_ = peak.mz_;
saristc/mzMLTrace
  data.resize(count);
saristc/mzMLTrace
  for (int i =0; i < count; i++) {
saristc/mzMLTrace
    data[i].mz_ = mzs[i];
saristc/mzMLTrace
    data[i].rt_ = rts[i];
saristc/mzMLTrace
    data[i].intn_ = intns[i];
saristc/mzMLTrace
  }
saristc/mzMLTrace
}
saristc/mzMLTrace

saristc/mzMLTrace
void specToChrom::get_chrom(int &ind, py::list &mz_s, py::list &intnss, 
saristc/mzMLTrace
                            py::list &rts){
saristc/mzMLTrace
  double hold;
saristc/mzMLTrace
  if (chroms_[ind].l_start_ == nullptr && chroms_[ind].h_start_ == nullptr) {
saristc/mzMLTrace
    for (int i = 0; i < chroms_[ind].chrom_len_; i++) {
saristc/mzMLTrace
      mz_s[i] = chroms_[ind].s_start_[i].mz_;
saristc/mzMLTrace
      intnss[i] = chroms_[ind].s_start_[i].intensity_;
saristc/mzMLTrace
      rts[i] = spectra_[chroms_[ind].s_start_[i].spec_id_].get_rt();
saristc/mzMLTrace
    }
saristc/mzMLTrace
  } else {
saristc/mzMLTrace
    std::fill(mz_quai_.begin(), mz_quai_.end() , 0.0);
saristc/mzMLTrace
    std::fill(in_quai_.begin(), in_quai_.end() , 0.0);
saristc/mzMLTrace
    std::fill(rt_quai_.begin(), rt_quai_.end() , 0.0);
saristc/mzMLTrace
    int it = 0;
saristc/mzMLTrace
    int count = 0;
saristc/mzMLTrace
    while ( count < chroms_[ind].l_count_ ) {
saristc/mzMLTrace
      if ( chroms_[ind].l_start_[it].chrom_id_ == ind ) {
saristc/mzMLTrace
        mz_quai_[ chroms_[ind].l_start_[it].spec_id_ - chroms_[ind].start_id_ ]
saristc/mzMLTrace
        += chroms_[ind].l_start_[it].mz_;
saristc/mzMLTrace
        in_quai_[ chroms_[ind].l_start_[it].spec_id_ - chroms_[ind].start_id_ ]
saristc/mzMLTrace
        += chroms_[ind].l_start_[it].intensity_;
saristc/mzMLTrace
        rt_quai_[ chroms_[ind].l_start_[it].spec_id_ - chroms_[ind].start_id_ ]
saristc/mzMLTrace
        = spectra_[chroms_[ind].l_start_[it].spec_id_].get_rt();
saristc/mzMLTrace
        count++;
saristc/mzMLTrace
      }
saristc/mzMLTrace
      it++;
saristc/mzMLTrace
    }
saristc/mzMLTrace
    it = 0;
saristc/mzMLTrace
    count = 0;
saristc/mzMLTrace
    while (count < chroms_[ind].s_count_) {
saristc/mzMLTrace
      if ( chroms_[ind].s_start_[it].chrom_id_ == ind ) {
saristc/mzMLTrace
        mz_quai_[ chroms_[ind].s_start_[it].spec_id_ - chroms_[ind].start_id_ ]
saristc/mzMLTrace
        += chroms_[ind].s_start_[it].mz_;
saristc/mzMLTrace
        in_quai_[ chroms_[ind].s_start_[it].spec_id_ - chroms_[ind].start_id_ ]
saristc/mzMLTrace
        += chroms_[ind].s_start_[it].intensity_;
saristc/mzMLTrace
        rt_quai_[ chroms_[ind].s_start_[it].spec_id_ - chroms_[ind].start_id_ ]
saristc/mzMLTrace
        = spectra_[chroms_[ind].s_start_[it].spec_id_].get_rt();
saristc/mzMLTrace
        count++;
saristc/mzMLTrace
      }
saristc/mzMLTrace
      it++;
saristc/mzMLTrace
    }
saristc/mzMLTrace
    it = 0;
saristc/mzMLTrace
    count = 0;
saristc/mzMLTrace
    while (count < chroms_[ind].h_count_) {
saristc/mzMLTrace
      if ( chroms_[ind].h_start_[it].chrom_id_ == ind ) {
saristc/mzMLTrace
        mz_quai_[ chroms_[ind].h_start_[it].spec_id_ - chroms_[ind].start_id_ ]
saristc/mzMLTrace
        += chroms_[ind].h_start_[it].mz_;
saristc/mzMLTrace
        in_quai_[ chroms_[ind].h_start_[it].spec_id_ - chroms_[ind].start_id_ ]
saristc/mzMLTrace
        += chroms_[ind].h_start_[it].intensity_;
saristc/mzMLTrace
        rt_quai_[ chroms_[ind].h_start_[it].spec_id_ - chroms_[ind].start_id_ ]
saristc/mzMLTrace
        = spectra_[chroms_[ind].h_start_[it].spec_id_].get_rt();
saristc/mzMLTrace
        count++;
saristc/mzMLTrace
      }
saristc/mzMLTrace
      it++;
saristc/mzMLTrace
    }
saristc/mzMLTrace
    for (int i = 0; i < chroms_[ind].chrom_len_; i++) {
saristc/mzMLTrace
      mz_s[i] = mz_quai_[i];
saristc/mzMLTrace
      intnss[i] = in_quai_[i];
saristc/mzMLTrace
      rts[i] = rt_quai_[i];
saristc/mzMLTrace
    }
saristc/mzMLTrace
  }
saristc/mzMLTrace
} // specToChrom::get_chrom
saristc/mzMLTrace

saristc/mzMLTrace
} // namespace asaristc
saristc/mzMLTrace

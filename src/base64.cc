#include <cstring>
saristc/mzMLTrace
#include <vector>
saristc/mzMLTrace

saristc/mzMLTrace
#include "base64.h"
saristc/mzMLTrace

saristc/mzMLTrace
namespace asaristc {
saristc/mzMLTrace

saristc/mzMLTrace
b64_decoder::b64_decoder() { form_table(); }
saristc/mzMLTrace

saristc/mzMLTrace
void b64_decoder::form_table() {
saristc/mzMLTrace
  if (tables_initiated_) {
saristc/mzMLTrace
    return;
saristc/mzMLTrace
  }
saristc/mzMLTrace

saristc/mzMLTrace
  int i, j, k;
saristc/mzMLTrace

saristc/mzMLTrace
  int prob_int = 1;
saristc/mzMLTrace
  is_little_endian_ = (*((char *)&prob_int) > 0 ? true : false);
saristc/mzMLTrace

saristc/mzMLTrace
  back_lookup_1_.resize(32767);
saristc/mzMLTrace
  back_lookup_2_.resize(32767);
saristc/mzMLTrace
  back_lookup_3_.resize(32767);
saristc/mzMLTrace

saristc/mzMLTrace
  main_lookup_.resize(16777214);
saristc/mzMLTrace

saristc/mzMLTrace
  for (i = '+'; i <= 'z'; i++) {
saristc/mzMLTrace
    for (j = '+'; j <= 'z'; j++) {
saristc/mzMLTrace
      int dest_adr = (i << 8) | j; /* could be +, but this is commutative */
saristc/mzMLTrace
      back_lookup_1_[dest_adr] = (base_table_[i] << 2) | (base_table_[j] >> 4);
saristc/mzMLTrace
      back_lookup_2_[dest_adr] = (base_table_[i] << 4) | (base_table_[j] >> 2);
saristc/mzMLTrace
      back_lookup_3_[dest_adr] = (base_table_[i] << 6) | (base_table_[j]);
saristc/mzMLTrace
    }
saristc/mzMLTrace
  }
saristc/mzMLTrace

saristc/mzMLTrace
  if (is_little_endian_) {
saristc/mzMLTrace
    for (i = '+'; i <= 'z'; i++) {
saristc/mzMLTrace
      for (j = '+'; j <= 'z'; j++) {
saristc/mzMLTrace
        for (k = '+'; k <= 'z'; k++) {
saristc/mzMLTrace
          int mailman = 0;
saristc/mzMLTrace
          char *c = (char *)&mailman;
saristc/mzMLTrace
          *c++ = i;
saristc/mzMLTrace
          *c++ = j;
saristc/mzMLTrace
          *c = k;
saristc/mzMLTrace
          mailman *= 2;
saristc/mzMLTrace
          main_lookup_[mailman++] = back_lookup_1_[(i << 8) | j];
saristc/mzMLTrace
          main_lookup_[mailman++] = back_lookup_2_[(j << 8) | k];
saristc/mzMLTrace
        }
saristc/mzMLTrace
      }
saristc/mzMLTrace
    }
saristc/mzMLTrace
  }
saristc/mzMLTrace
}
saristc/mzMLTrace

saristc/mzMLTrace
void b64_decoder::decode_base64(char *bin_data, const char *base64_data,
saristc/mzMLTrace
                                long unsigned int &bytes_to_decode) {
saristc/mzMLTrace
  unsigned char *bin_dest = (unsigned char *)bin_data;
saristc/mzMLTrace
  const char *char_src = base64_data;
saristc/mzMLTrace
  int bytes_left = static_cast<int>(bytes_to_decode);
saristc/mzMLTrace
  unsigned short int f;
saristc/mzMLTrace
  int data_add = 0;
saristc/mzMLTrace
  while (bytes_left >= 3) {
saristc/mzMLTrace
    memcpy((char *)&data_add, char_src, 3);
saristc/mzMLTrace
    memcpy(bin_dest, &main_lookup_[data_add * 2], 2);
saristc/mzMLTrace
    *(bin_dest + 2) =
saristc/mzMLTrace
        back_lookup_3_[((*(char_src + 2)) << 8) | (*(char_src + 3))];
saristc/mzMLTrace
    char_src += 4;
saristc/mzMLTrace
    bin_dest += 3;
saristc/mzMLTrace
    bytes_left -= 3;
saristc/mzMLTrace
  }
saristc/mzMLTrace

saristc/mzMLTrace
  while (bytes_left--) {
saristc/mzMLTrace
    *bin_dest++ = back_lookup_1_[f = (((*char_src) << 8) | (*(char_src + 1)))];
saristc/mzMLTrace
    if (bytes_left--) {
saristc/mzMLTrace
      *bin_dest++ = back_lookup_2_[f = ((f << 8) | (*(char_src + 2)))];
saristc/mzMLTrace
      if (bytes_left--) {
saristc/mzMLTrace
        *bin_dest++ = back_lookup_3_[(short int)((f << 8) | (*(char_src + 3)))];
saristc/mzMLTrace
      } else {
saristc/mzMLTrace
        break;
saristc/mzMLTrace
      }
saristc/mzMLTrace
    } else {
saristc/mzMLTrace
      break;
saristc/mzMLTrace
    }
saristc/mzMLTrace
    char_src += 4;
saristc/mzMLTrace
  }
saristc/mzMLTrace
}
saristc/mzMLTrace

saristc/mzMLTrace
} // namespace asaristc
saristc/mzMLTrace

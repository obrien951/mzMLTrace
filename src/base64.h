#ifndef _ASARI_BASE64_
saristc/mzMLTrace
#define _ASARI_BASE64_
saristc/mzMLTrace
#include <vector>
saristc/mzMLTrace

saristc/mzMLTrace
namespace asaristc {
saristc/mzMLTrace

saristc/mzMLTrace
class b64_decoder {
saristc/mzMLTrace
public:
saristc/mzMLTrace
  b64_decoder();
saristc/mzMLTrace

saristc/mzMLTrace
  //  void decode( unsigned char * &inp, unsigned char * &outp);
saristc/mzMLTrace
  void decode_base64(char *bin_data, const char *base64_data,
saristc/mzMLTrace
                     long unsigned int &bytes_to_decode);
saristc/mzMLTrace

saristc/mzMLTrace
protected:
saristc/mzMLTrace
  // Checks to see if the processor is bigendian or littlendian then forms
saristc/mzMLTrace
  //   the appropriate look-up table. The table isn't hard-coded to lower the
saristc/mzMLTrace
  //   memory footprint
saristc/mzMLTrace
  void form_table();
saristc/mzMLTrace

saristc/mzMLTrace
  bool tables_initiated_ = false;
saristc/mzMLTrace
  bool is_little_endian_;
saristc/mzMLTrace

saristc/mzMLTrace
  /* Look up tables for base64 to binary conversion */
saristc/mzMLTrace
  /* Copy values out of this table 3 bits at a time.
saristc/mzMLTrace
   * Decodes the bulk of the data */
saristc/mzMLTrace
  std::vector<unsigned char> main_lookup_;
saristc/mzMLTrace

saristc/mzMLTrace
  /* Dealing with padded data. Are used to decode the first second and third
saristc/mzMLTrace
   * dangling characters respectively */
saristc/mzMLTrace
  std::vector<unsigned char> back_lookup_1_;
saristc/mzMLTrace
  std::vector<unsigned char> back_lookup_2_;
saristc/mzMLTrace
  std::vector<unsigned char> back_lookup_3_;
saristc/mzMLTrace

saristc/mzMLTrace
  unsigned int base_table_[128] = {
saristc/mzMLTrace
      // basic base64 charset table
saristc/mzMLTrace
      0,  //  NUL
saristc/mzMLTrace
      0,  //  SOH
saristc/mzMLTrace
      0,  //  STX
saristc/mzMLTrace
      0,  //  ETX
saristc/mzMLTrace
      0,  //  EOT
saristc/mzMLTrace
      0,  //  ENQ
saristc/mzMLTrace
      0,  //  ACK
saristc/mzMLTrace
      0,  //  BEL
saristc/mzMLTrace
      0,  //   BS
saristc/mzMLTrace
      0,  //   HT
saristc/mzMLTrace
      0,  //   LF
saristc/mzMLTrace
      0,  //   VT
saristc/mzMLTrace
      0,  //   FF
saristc/mzMLTrace
      0,  //   CR
saristc/mzMLTrace
      0,  //   SO
saristc/mzMLTrace
      0,  //   SI
saristc/mzMLTrace
      0,  //  DLE
saristc/mzMLTrace
      0,  //  DC1
saristc/mzMLTrace
      0,  //  DC2
saristc/mzMLTrace
      0,  //  DC3
saristc/mzMLTrace
      0,  //  DC4
saristc/mzMLTrace
      0,  //  NAK
saristc/mzMLTrace
      0,  //  SYN
saristc/mzMLTrace
      0,  //  ETB
saristc/mzMLTrace
      0,  //  CAN
saristc/mzMLTrace
      0,  //   EM
saristc/mzMLTrace
      0,  //  SUB
saristc/mzMLTrace
      0,  //  ESC
saristc/mzMLTrace
      0,  //   FS
saristc/mzMLTrace
      0,  //   GS
saristc/mzMLTrace
      0,  //   RS
saristc/mzMLTrace
      0,  //   US
saristc/mzMLTrace
      0,  //   SP
saristc/mzMLTrace
      0,  //    !
saristc/mzMLTrace
      0,  //    "
saristc/mzMLTrace
      0,  //    #
saristc/mzMLTrace
      0,  //    $
saristc/mzMLTrace
      0,  //    %
saristc/mzMLTrace
      0,  //    &
saristc/mzMLTrace
      0,  //    '
saristc/mzMLTrace
      0,  //    (
saristc/mzMLTrace
      0,  //    )
saristc/mzMLTrace
      0,  //    *
saristc/mzMLTrace
      62, //    +
saristc/mzMLTrace
      0,  //    ,
saristc/mzMLTrace
      0,  //    -
saristc/mzMLTrace
      0,  //    .
saristc/mzMLTrace
      63, //    /
saristc/mzMLTrace
      52, //    0,
saristc/mzMLTrace
      53, //    1
saristc/mzMLTrace
      54, //    2
saristc/mzMLTrace
      55, //    3
saristc/mzMLTrace
      56, //    4
saristc/mzMLTrace
      57, //    5
saristc/mzMLTrace
      58, //    6
saristc/mzMLTrace
      59, //    7
saristc/mzMLTrace
      60, //    8
saristc/mzMLTrace
      61, //    9
saristc/mzMLTrace
      0,  //    :
saristc/mzMLTrace
      0,  //    ;
saristc/mzMLTrace
      0,  //    <
saristc/mzMLTrace
      0,  //    =
saristc/mzMLTrace
      0,  //    >
saristc/mzMLTrace
      0,  //    ?
saristc/mzMLTrace
      0,  //    @
saristc/mzMLTrace
      0,  //    A
saristc/mzMLTrace
      1,  //    B
saristc/mzMLTrace
      2,  //    C
saristc/mzMLTrace
      3,  //    D
saristc/mzMLTrace
      4,  //    E
saristc/mzMLTrace
      5,  //    F
saristc/mzMLTrace
      6,  //    G
saristc/mzMLTrace
      7,  //    H
saristc/mzMLTrace
      8,  //    I
saristc/mzMLTrace
      9,  //    J
saristc/mzMLTrace
      10, //    K
saristc/mzMLTrace
      11, //    L
saristc/mzMLTrace
      12, //    M
saristc/mzMLTrace
      13, //    N
saristc/mzMLTrace
      14, //    O
saristc/mzMLTrace
      15, //    P
saristc/mzMLTrace
      16, //    Q
saristc/mzMLTrace
      17, //    R
saristc/mzMLTrace
      18, //    S
saristc/mzMLTrace
      19, //    T
saristc/mzMLTrace
      20, //    U
saristc/mzMLTrace
      21, //    V
saristc/mzMLTrace
      22, //    W
saristc/mzMLTrace
      23, //    X
saristc/mzMLTrace
      24, //    Y
saristc/mzMLTrace
      25, //    Z
saristc/mzMLTrace
      0,  //    [
saristc/mzMLTrace
      0,  //    '\'
saristc/mzMLTrace
      0,  //    ]
saristc/mzMLTrace
      0,  //    ^
saristc/mzMLTrace
      0,  //    _
saristc/mzMLTrace
      0,  //    `
saristc/mzMLTrace
      26, //    a
saristc/mzMLTrace
      27, //    b
saristc/mzMLTrace
      28, //    c
saristc/mzMLTrace
      29, //    d
saristc/mzMLTrace
      30, //    e
saristc/mzMLTrace
      31, //    f
saristc/mzMLTrace
      32, //    g
saristc/mzMLTrace
      33, //    h
saristc/mzMLTrace
      34, //    i
saristc/mzMLTrace
      35, //    j
saristc/mzMLTrace
      36, //    k
saristc/mzMLTrace
      37, //    l
saristc/mzMLTrace
      38, //    m
saristc/mzMLTrace
      39, //    n
saristc/mzMLTrace
      40, //    o
saristc/mzMLTrace
      41, //    p
saristc/mzMLTrace
      42, //    q
saristc/mzMLTrace
      43, //    r
saristc/mzMLTrace
      44, //    s
saristc/mzMLTrace
      45, //    t
saristc/mzMLTrace
      46, //    u
saristc/mzMLTrace
      47, //    v
saristc/mzMLTrace
      48, //    w
saristc/mzMLTrace
      49, //    x
saristc/mzMLTrace
      50, //    y
saristc/mzMLTrace
      51, //    z
saristc/mzMLTrace
      0,  //    {
saristc/mzMLTrace
      0,  //    |
saristc/mzMLTrace
      0,  //    }
saristc/mzMLTrace
      0,  //    ~
saristc/mzMLTrace
      0   //  DEL
saristc/mzMLTrace
  };
saristc/mzMLTrace
};
saristc/mzMLTrace

saristc/mzMLTrace
} // namespace asaristc
saristc/mzMLTrace

saristc/mzMLTrace
#else
saristc/mzMLTrace
#endif
saristc/mzMLTrace

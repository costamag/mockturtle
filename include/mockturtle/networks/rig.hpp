/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYrigHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file rig.hpp
  \brief Representation Independent Graph

  This network assumes that buffers, inverter and splitters are cost free.
  Everything you declare apart for these has a cost.
  The network is structurally hashed for gates of the same type.
  gates of different type are not hashed together even if related by negation.
  create_and(x1,x2) != !create_nand(x1,x2).
  but naturally reate_and(x1,x2) == !create_and(x1,x2).
  Any overwriting is a representation-dependent assumption.

  \author Andrea Costamagna
*/

#pragma once

#include "../utils/truth_table_cache.hpp"
#include "../utils/algorithm.hpp"
#include "detail/foreach.hpp"
#include "../utils/tech_library.hpp"
#include "../algorithms/node_resynthesis/xag_npn.hpp"
#include "aig.hpp"
#include <kitty/constructors.hpp>
#include <kitty/operations.hpp>
#include "storage.hpp"
#include <kitty/kitty.hpp>
#include <algorithm>
#include <memory>
#include <list>

namespace mockturtle
{

namespace rils
{

/*! \brief identifier of the node category */
enum i_func_t : uint32_t
{
  i_CONST = 0u,
  i_PI = 1u,
  i_BUF = 2u
};

/*! \brief literals of the precomputed truth-tables in the truth-table cache.
  *
  * The default convention for literals is assumed.  That is an index \f$i\f$
  * (starting) from \f$0\f$ has positive literal \f$2i\f$ and negative literal
  * \f$2i + 1\f$.
  * 
  */
enum e_func_t : uint32_t
{
  e_CONST = 0u,
  e_PI = 1u,
  e_BUF = 2u,
  e_AND = 4u,
  e_OR = 6u,
  e_LT = 8u,
  e_GT = 10u,
  e_XOR = 12u
};

#pragma region utils
  template<typename PTR_T>
  struct signal_t
  {
    signal_t() = default;

    signal_t( uint64_t index, uint64_t complement )
        : complement( complement ), index( index )
    {
    }

    signal_t( uint32_t index )
        : complement( 0 ), index( index )
    {
    }

    signal_t( uint64_t index, uint64_t complement, uint64_t output )
        : complement( complement ), index( index )
    {
    }

    explicit signal_t( uint64_t data )
        : data( data )
    {
    }

    signal_t( PTR_T const& p )
        : complement( p.weight & 1 ), index( p.index )
    {
    }

    union
    {
      struct
      {
        uint64_t complement : 1;
        uint64_t index : 63;
      };
      uint64_t data;
    };

    signal_t operator!() const
    {
      return signal_t( data ^ 1 );
    }

    signal_t operator+() const
    {
      return { index, 0 };
    }

    signal_t operator-() const
    {
      return { index, 1 };
    }

    signal_t operator^( bool complement ) const
    {
      return signal_t( data ^ ( complement ? 1 : 0 ) );
    }

    bool operator==( signal_t const& other ) const
    {
      return data == other.data;
    }

    bool operator!=( signal_t const& other ) const
    {
      return data != other.data;
    }

    bool operator<( signal_t const& other ) const
    {
      return data < other.data;
    }

    operator PTR_T() const
    {
      return { index, complement };
    }

    operator uint64_t() const
    {
      return data;
    }

  #if __cplusplus > 201703L
    bool operator==( PTR_T const& other ) const
    {
      return data == other.data;
    }
  #endif
  };

  /*! \brief rig luts storage
  */
  struct e_data_t
  {
    truth_table_cache<kitty::dynamic_truth_table> cache;
  };

  struct e_gate_t
  {
    using pointer_type = node_pointer<1>;
    using twin_pointer_type = node_pointer<1>;

    std::vector<pointer_type> children;

    /*! \brief number of fanouts */
    uint32_t nfos{0};
    /*! \brief id of the functionality stored in the tt-cache */
    uint32_t func{0};
    /*! \brief application specific bits */
    uint32_t bits{0};
    /*! \brief Twin internal signal */
    signal_t<twin_pointer_type> twin {0,0};

    bool operator==( e_gate_t const& other ) const
    {
      return (func == other.func) && (children == other.children);
    }
  };

  struct i_gate_t
  {
    using pointer_type = node_pointer<1>;
    using twin_pointer_type = node_pointer<1>;

    std::array<pointer_type, 2> children;
    /*! \brief number of fanouts */
    uint32_t nfos{0};
    /*! \brief node type const:0 pi:1 */
    uint32_t func{0};
    /*! \brief application specific bits */
    uint32_t bits{0};
    /*! Twin external signal: W might not exit! */
    signal_t<twin_pointer_type> twin {0,0};

    bool operator==( i_gate_t const& other ) const
    {
      return children == other.children;
    }

  };

  /*! \brief Hash function for AIGs (from ABC) */
  struct i_hash_t
  {
    uint64_t operator()( i_gate_t const& n ) const
    {
      uint64_t seed = -2011;
      seed += n.children[0].index * 7937;
      seed += n.children[1].index * 2971;
      seed += n.children[0].weight * 911;
      seed += n.children[1].weight * 353;
      return seed;
    }
  };
  using e_storage_t = smart_storage<e_gate_t, e_data_t>;

  using i_data_t = empty_storage_data;
  struct i_hash_t;
  using i_storage_t = smart_storage<i_gate_t, i_data_t, i_hash_t>;
#pragma endregion utils

class rig_network
{

  #pragma region types
  public:
    using e_node_t = uint64_t;
    using e_signal_t = signal_t<e_gate_t::pointer_type>;

    using i_node_t = uint64_t;
    using i_signal_t = signal_t<i_gate_t::pointer_type>;

    static constexpr auto min_fanin_size = 1;
    static constexpr auto max_fanin_size = 32;
    using base_t = rig_network;
    using node = e_node_t;
    using signal = e_signal_t;
  #pragma endregion types

  #pragma region constructor
  public:
    rig_network();

  protected:
    inline void _init();
  #pragma endregion constructors

  #pragma region linking
  private:
    e_signal_t get_e_signal( i_signal_t const& );
    i_signal_t get_i_signal( e_signal_t const& );
  #pragma endregion linking

  #pragma region Primary I / O and constants
  public:
    signal get_constant( bool );
    signal create_pi();
    uint32_t create_po( signal const& );
    bool is_combinational() const;
    bool is_constant( node const& ) const;
    bool is_ci( node const& ) const;
    bool is_pi( node const& ) const;
  #pragma endregion Primary I / O and constants

  #pragma region nodes and signals
  public:
    node get_node( signal const& ) const;
    signal make_signal( node const& ) const;
    bool is_complemented( signal const& ) const;
    uint32_t node_to_index( node const& ) const;
    node index_to_node( uint32_t ) const;
    node pi_at( uint32_t ) const;
    node ci_at( uint32_t ) const;
    signal po_at( uint32_t ) const;
    signal co_at( uint32_t ) const;
    uint32_t ci_index( node const& ) const;
    uint32_t pi_index( node const& ) const;
    //uint32_t co_index( signal const& ) const;
    //uint32_t po_index( signal const& ) const;
  #pragma endregion nodes and signals

  #pragma region node and signal iterators
    template<typename Fn> void foreach_node( Fn&& fn ) const;
    template<typename Fn> void foreach_ci( Fn&& fn ) const;
    template<typename Fn> void foreach_co( Fn&& fn ) const;
    template<typename Fn> void foreach_pi( Fn&& fn ) const;
    template<typename Fn> void foreach_po( Fn&& fn ) const;
    template<typename Fn> void foreach_gate( Fn&& fn ) const;
    template<typename Fn> void foreach_fanin( node const& n, Fn&& fn ) const;
  #pragma endregion node and signal iterators

  #pragma region unary functions
    signal create_buf( signal const& );
    signal create_not( signal const& );
    bool is_buf( node const& );
    bool is_not( node const& );
  #pragma endregion unary functions

  #pragma region binary functions
    i_signal_t i_create_and( i_signal_t, i_signal_t );
    i_signal_t i_create_xor( i_signal_t, i_signal_t );
    signal create_and( signal, signal );
    signal create_nand( signal, signal );
    signal create_or( signal, signal );
    signal create_nor( signal, signal );
    signal create_lt( signal, signal );
    signal create_ge( signal, signal );
    signal create_gt( signal, signal );
    signal create_le( signal, signal );
    signal create_xor( signal, signal );
    signal create_xnor( signal, signal );
    bool is_and( node const& );
    bool is_nand( node const& );
    bool is_or( node const& );
    bool is_nor( node const& );
    bool is_lt( node const& );
    bool is_ge( node const& );
    bool is_gt( node const& );
    bool is_le( node const& );
    bool is_xor( node const& );
    bool is_xnor( node const& );
  #pragma endregion binary functions

  #pragma region arbitrary function
    signal _create_node( std::vector<signal> const&, uint32_t );
    std::tuple<rig_network::signal, bool> _create_known_node( std::vector<signal> const&, uint32_t );
    bool is_function( node const& ) const;

    i_signal_t synthesize_twin( std::array<i_signal_t, 4u> &, uint32_t );
    i_signal_t synthesize_twin_rec( std::array<i_signal_t, 4u> &, kitty::dynamic_truth_table const& );
    i_signal_t boolean_matching( std::array<i_signal_t, 4u> &, kitty::dynamic_truth_table const& );
  #pragma endregion arbitrary function

  #pragma region structural properties
  public:
    size_t size() const;
    size_t num_cis() const;
    size_t num_cos() const;
    size_t num_pis() const;
    size_t num_pos() const;
    size_t num_gates() const;
    size_t fanin_size( node const& n ) const;
    size_t fanout_size( node const& n ) const;
    size_t incr_fanout_size( node const& n ) const;
    size_t decr_fanout_size( node const& n ) const;
  #pragma endregion structural properties

  #pragma region functional properties
  public:
    kitty::dynamic_truth_table node_function( const node& ) const;
  #pragma endregion functional properties

public:
  std::shared_ptr<e_storage_t> _e_storage;
  std::shared_ptr<i_storage_t> _i_storage;
  /* complete AIG database */
  inline static const uint32_t subgraphs_aig [] = { 52223492, 4, 7, 5,
    6, 11, 13, 8, 14, 9, 15, 17, 19, 3, 5, 6, 9, 7, 8, 25, 27, 23, 29,
    22, 28, 31, 33, 3, 4, 2, 6, 37, 39, 8, 40, 9, 41, 43, 45, 2, 5, 2,
    7, 4, 51, 49, 53, 8, 55, 9, 54, 57, 59, 2, 26, 25, 63, 2, 63, 5,
    67, 65, 68, 64, 69, 71, 73, 2, 4, 6, 23, 77, 79, 8, 80, 9, 81, 83,
    85, 5, 39, 8, 88, 2, 8, 6, 39, 93, 95, 89, 96, 91, 99, 4, 8, 2,
    103, 13, 104, 13, 27, 3, 109, 107, 111, 8, 76, 6, 22, 77, 117, 27,
    118, 115, 121, 7, 9, 22, 125, 6, 8, 23, 129, 76, 125, 130, 133,
    127, 135, 5, 8, 3, 139, 93, 141, 4, 141, 7, 145, 142, 146, 143,
    147, 149, 151, 5, 7, 8, 155, 9, 154, 157, 159, 2, 161, 4, 6, 161,
    165, 3, 167, 163, 169, 2, 9, 4, 173, 3, 8, 7, 176, 174, 179, 39,
    179, 5, 183, 181, 185, 4, 9, 155, 189, 2, 190, 165, 190, 3, 195,
    193, 197, 23, 77, 6, 77, 8, 203, 200, 205, 201, 204, 207, 209, 7,
    23, 77, 212, 9, 76, 23, 217, 6, 219, 215, 221, 23, 25, 77, 129,
    224, 226, 225, 227, 229, 231, 3, 9, 7, 103, 235, 236, 234, 237,
    239, 241, 23, 115, 7, 245, 8, 23, 77, 249, 6, 251, 247, 253, 23,
    125, 227, 257, 226, 256, 259, 261, 6, 200, 7, 201, 265, 267, 9,
    269, 49, 177, 7, 273, 25, 275, 7, 173, 4, 177, 278, 281, 25, 283,
    7, 77, 8, 287, 125, 289, 9, 165, 155, 293, 3, 295, 8, 164, 2, 294,
    299, 301, 297, 302, 4, 50, 5, 234, 307, 309, 129, 310, 11, 23,
    177, 315, 129, 317, 37, 49, 6, 320, 7, 321, 323, 325, 8, 326, 9,
    327, 329, 331, 9, 77, 7, 334, 23, 337, 129, 339, 9, 286, 26, 77,
    9, 23, 6, 346, 345, 349, 4, 234, 7, 139, 2, 354, 353, 357, 129,
    358, 77, 256, 129, 362, 2, 124, 4, 366, 3, 125, 4, 371, 367, 373,
    129, 374, 369, 377, 203, 249, 6, 76, 343, 383, 289, 384, 8, 325,
    9, 324, 389, 391, 8, 11, 51, 394, 51, 235, 4, 399, 397, 401, 367,
    371, 4, 405, 5, 404, 407, 409, 129, 411, 3, 6, 249, 415, 5, 176,
    4, 129, 2, 420, 419, 423, 6, 424, 7, 425, 427, 429, 5, 9, 2, 432,
    157, 435, 3, 7, 7, 439, 2, 441, 4, 443, 439, 444, 8, 445, 441,
    448, 447, 451, 9, 383, 24, 37, 325, 457, 9, 37, 11, 461, 50, 463,
    51, 462, 465, 467, 88, 439, 173, 439, 4, 473, 471, 475, 50, 139,
    139, 189, 3, 481, 479, 483, 129, 484, 155, 235, 293, 488, 292,
    489, 491, 493, 6, 321, 9, 51, 320, 499, 497, 501, 5, 172, 129,
    505, 4, 278, 506, 509, 7, 200, 24, 201, 513, 515, 7, 22, 347, 519,
    337, 521, 5, 399, 398, 420, 525, 527, 22, 286, 287, 346, 531, 533,
    3, 129, 4, 536, 2, 128, 4, 125, 541, 543, 537, 544, 539, 547, 8,
    77, 212, 551, 9, 201, 213, 555, 553, 557, 154, 370, 9, 371, 155,
    562, 561, 565, 24, 320, 26, 321, 569, 571, 8, 286, 9, 287, 575,
    577, 2, 542, 3, 543, 7, 583, 4, 543, 9, 587, 585, 589, 581, 591,
    25, 325, 125, 129, 9, 212, 8, 213, 77, 601, 598, 602, 599, 603,
    605, 607, 7, 93, 4, 235, 611, 612, 610, 613, 615, 617, 129, 618,
    5, 370, 7, 623, 8, 625, 373, 627, 3, 236, 129, 631, 77, 632, 125,
    201, 129, 637, 3, 154, 8, 641, 155, 165, 234, 645, 643, 647, 5,
    124, 3, 651, 542, 653, 543, 652, 129, 657, 655, 658, 4, 370, 5,
    367, 371, 664, 663, 667, 129, 669, 173, 641, 3, 124, 5, 674, 581,
    677, 129, 678, 9, 49, 155, 682, 154, 683, 685, 687, 248, 286, 22,
    287, 76, 249, 693, 695, 691, 696, 9, 22, 265, 701, 115, 702, 334,
    415, 157, 707, 6, 460, 4, 37, 7, 713, 8, 714, 711, 717, 2, 165, 9,
    721, 155, 722, 154, 723, 725, 727, 8, 154, 9, 155, 731, 733, 5,
    129, 2, 736, 173, 737, 739, 741, 165, 172, 157, 745, 13, 103, 51,
    188, 155, 751, 9, 50, 748, 755, 4, 24, 155, 759, 7, 138, 759, 763,
    2, 645, 9, 767, 731, 769, 8, 720, 6, 721, 9, 774, 773, 777, 4,
    778, 5, 779, 781, 783, 12, 173, 6, 172, 9, 789, 4, 791, 787, 793,
    2, 237, 4, 796, 129, 799, 631, 800, 125, 219, 7, 218, 77, 806,
    805, 809, 6, 235, 9, 813, 8, 812, 5, 817, 815, 819, 309, 821, 92,
    154, 3, 164, 93, 827, 155, 828, 825, 831, 188, 789, 279, 789, 5,
    837, 835, 839, 9, 644, 9, 39, 5, 844, 737, 845, 847, 849, 3, 138,
    25, 51, 2, 139, 5, 856, 855, 859, 853, 861, 77, 433, 7, 865, 9,
    865, 6, 869, 867, 871, 5, 439, 175, 875, 234, 293, 235, 295, 879,
    881, 188, 836, 839, 885, 27, 77, 25, 201, 889, 891, 888, 890, 893,
    895, 235, 237, 7, 234, 5, 900, 103, 903, 813, 904, 5, 536, 125,
    909, 6, 234, 9, 235, 4, 915, 7, 917, 913, 919, 8, 439, 237, 922,
    236, 923, 925, 927, 4, 172, 6, 930, 7, 852, 933, 935, 53, 682,
    157, 165, 6, 37, 8, 321, 943, 945, 129, 947, 5, 23, 6, 951, 249,
    953, 2, 154, 9, 957, 165, 958, 731, 961, 9, 323, 154, 234, 489,
    967, 160, 165, 5, 415, 9, 972, 7, 972, 8, 977, 975, 979, 8, 518,
    9, 519, 983, 985, 22, 26, 347, 989, 8, 320, 497, 993, 643, 827, 3,
    165, 2, 155, 999, 1001, 157, 1003, 3, 433, 154, 1007, 6, 1006, 9,
    1011, 155, 1013, 1009, 1015, 3, 11, 8, 1019, 27, 1018, 1021, 1023,
    8, 640, 173, 1027, 7, 235, 5, 1030, 4, 1031, 9, 1034, 1033, 1037,
    4, 537, 6, 1041, 9, 1043, 909, 1045, 6, 177, 5, 1049, 179, 1050,
    173, 179, 4, 1055, 1053, 1057, 9, 999, 5, 998, 7, 1062, 1061,
    1065, 129, 675, 581, 1068, 3, 13, 8, 1072, 5, 1074, 2, 12, 9,
    1073, 1079, 1080, 1077, 1083, 7, 857, 4, 857, 6, 1089, 9, 1090,
    1087, 1093, 324, 461, 325, 460, 1097, 1099, 93, 125, 415, 1102, 5,
    1105, 543, 1107, 5, 93, 7, 1110, 6, 1111, 2, 1114, 3, 1115, 1117,
    1119, 9, 1121, 1113, 1123, 8, 998, 2, 164, 9, 1128, 155, 1131,
    1126, 1132, 1127, 1133, 1135, 1137, 6, 335, 249, 701, 1141, 1142,
    77, 701, 7, 1147, 6, 1146, 249, 1151, 1149, 1152, 23, 579, 22,
    578, 1157, 1159, 77, 155, 93, 1163, 4, 415, 1031, 1167, 4, 1030,
    93, 1171, 1169, 1172, 9, 415, 155, 1176, 641, 1179, 287, 383, 249,
    1183, 125, 488, 124, 489, 1187, 1189, 5, 414, 2, 10, 1193, 1195,
    8, 1197, 9, 1196, 1199, 1201, 155, 172, 154, 176, 1205, 1207, 9,
    766, 730, 767, 1211, 1213, 875, 1007, 129, 1217, 455, 519, 23,
    165, 6, 173, 1223, 1225, 383, 641, 9, 1229, 643, 1231, 6, 201,
    249, 1235, 93, 415, 6, 415, 5, 1241, 1238, 1243, 8, 1242, 1239,
    1246, 1245, 1249, 155, 723, 9, 973, 53, 1254, 454, 519, 2, 293, 3,
    292, 157, 1263, 1261, 1264, 23, 1001, 9, 1268, 7, 1268, 8, 1273,
    1271, 1275, 5, 371, 37, 1279, 129, 1281, 165, 643, 49, 103, 7,
    1286, 6, 1287, 9, 1290, 1289, 1293, 7, 335, 265, 1297, 249, 1299,
    1050, 1054, 1057, 1303, 129, 200, 155, 1129, 93, 1309, 22, 124, 3,
    102, 49, 1315, 7, 1317, 25, 1319, 9, 645, 643, 1323, 9, 164, 154,
    173, 1327, 1329, 9, 640, 1284, 1333, 286, 347, 533, 1337, 27, 235,
    10, 1341, 11, 1340, 1343, 1345, 2, 759, 9, 1348, 3, 761, 1351,
    1353, 461, 714, 460, 715, 1357, 1359, 9, 537, 7, 1362, 909, 1365,
    23, 423, 596, 1368, 597, 1369, 1371, 1373, 663, 1279, 129, 1377,
    641, 733, 165, 1381, 3, 189, 7, 1385, 4, 1386, 9, 1385, 6, 1391,
    1389, 1393, 9, 79, 155, 1061, 1129, 1399, 226, 257, 1313, 1402, 8,
    22, 286, 1407, 221, 1409, 9, 1129, 155, 998, 1412, 1415, 4, 439,
    173, 1418, 439, 1049, 5, 1423, 1421, 1425, 160, 1129, 7, 37, 9,
    36, 1431, 1433, 49, 1435, 129, 1437, 5, 51, 1177, 1441, 1176,
    1440, 129, 1445, 1443, 1446, 9, 213, 212, 550, 1451, 1453, 157,
    1129, 3, 158, 1456, 1459, 3, 155, 2, 733, 1463, 1465, 39, 460,
    235, 875, 234, 874, 1471, 1473, 103, 1475, 125, 908, 124, 909,
    1479, 1481, 675, 911, 2, 189, 51, 737, 1487, 1489, 4, 124, 909,
    1493, 675, 1495, 165, 1463, 9, 1499, 731, 1501, 6, 461, 3, 460, 5,
    1507, 1430, 1509, 1505, 1511, 2, 188, 23, 1515, 7, 1517, 6, 1516,
    9, 1520, 1519, 1523, 292, 957, 5, 438, 581, 1529, 129, 1530, 79,
    519, 9, 1534, 7, 347, 77, 1538, 349, 1541, 80, 984, 226, 519, 77,
    415, 9, 1549, 155, 1551, 79, 334, 8, 644, 165, 1557, 2, 1559,
    1323, 1557, 3, 1562, 1561, 1565, 129, 373, 23, 1568, 6, 433, 3,
    1573, 8, 1575, 5, 433, 1574, 1579, 1577, 1581, 9, 13, 2, 1585, 8,
    12, 5, 1589, 3, 1591, 1587, 1593, 8, 48, 1430, 1597, 1505, 1599,
    20, 35, 47, 60, 74, 87, 101, 113, 122, 137, 152, 170, 186, 198,
    210, 223, 233, 243, 254, 263, 270, 277, 285, 290, 304, 312, 318,
    333, 340, 342, 351, 360, 364, 379, 380, 386, 392, 402, 412, 416,
    430, 436, 452, 454, 459, 468, 477, 390, 486, 494, 502, 510, 517,
    522, 529, 535, 549, 558, 567, 573, 579, 592, 595, 596, 608, 620,
    628, 634, 638, 648, 660, 670, 673, 680, 689, 699, 704, 708, 719,
    729, 735, 742, 746, 748, 753, 756, 761, 765, 771, 784, 794, 802,
    811, 822, 833, 292, 841, 842, 850, 863, 872, 876, 883, 887, 897,
    899, 906, 911, 921, 928, 937, 938, 489, 940, 948, 954, 963, 964,
    968, 970, 980, 987, 9, 991, 994, 996, 1004, 1016, 1024, 1029,
    1039, 1047, 1059, 1067, 1070, 1085, 1095, 1101, 1108, 1125, 1139,
    1144, 1154, 1161, 1164, 1174, 1181, 1184, 1190, 1203, 1209, 1215,
    1218, 1221, 1226, 1232, 1236, 1251, 1253, 521, 1256, 1258, 1266,
    1276, 1282, 1284, 1295, 1300, 1305, 1306, 1310, 1312, 1321, 1324,
    1331, 1334, 1339, 1346, 1355, 1361, 1367, 1375, 1378, 1382, 1394,
    1396, 1400, 1404, 1411, 1416, 1426, 1428, 1438, 1448, 1455, 1460,
    1466, 1468, 1476, 1483, 1484, 1490, 1496, 1503, 124, 158, 1512, 8,
    1525, 1526, 1532, 1536, 1543, 1544, 1546, 1553, 1554, 1566, 1570,
    1582, 1594, 1600 };

};


#pragma region constructors
  /*! \brief Network constructor

  * Construct the network using the init function
  */
  rig_network::rig_network()
    : _e_storage( std::make_shared<e_storage_t>() ), _i_storage( std::make_shared<i_storage_t>() )
  {
    _init();
  }
  /*! \brief Network initializer

  * At initialization, the network must have allocated only one node for constant 0.
  * This method stores the truth table of the constant function and connects the constant 0
  * node of the externale and the internal representations.
  */
  inline void rig_network::_init()
  {
    /* already initialized */
    if ( _e_storage->nodes.size() > 1 )
      return;

    /* constant node : #0 in the cache */
    kitty::dynamic_truth_table tt_zero( 0 );
    _e_storage->data.cache.insert( tt_zero );
    _i_storage->nodes[0].func = i_func_t::i_CONST;

    static uint64_t _not = 0x1;
    kitty::dynamic_truth_table tt_not( 1 );
    kitty::create_from_words( tt_not, &_not, &_not + 1 );
    _e_storage->data.cache.insert( tt_not );

    static uint64_t _and = 0x8;
    kitty::dynamic_truth_table tt_and( 2 );
    kitty::create_from_words( tt_and, &_and, &_and + 1 );
    _e_storage->data.cache.insert( tt_and );

    static uint64_t _or = 0xe;
    kitty::dynamic_truth_table tt_or( 2 );
    kitty::create_from_words( tt_or, &_or, &_or + 1 );
    _e_storage->data.cache.insert( tt_or );

    static uint64_t _lt = 0x2;
    kitty::dynamic_truth_table tt_lt( 2 );
    kitty::create_from_words( tt_lt, &_lt, &_lt + 1 );
    _e_storage->data.cache.insert( tt_lt );

    static uint64_t _gt = 0x4;
    kitty::dynamic_truth_table tt_gt( 2 );
    kitty::create_from_words( tt_gt, &_gt, &_gt + 1 );
    _e_storage->data.cache.insert( tt_gt );

    static uint64_t _xor = 0x6;
    kitty::dynamic_truth_table tt_xor( 2 );
    kitty::create_from_words( tt_xor, &_xor, &_xor + 1 );
    _e_storage->data.cache.insert( tt_xor );

    static uint64_t _maj = 0xe8;
    kitty::dynamic_truth_table tt_maj( 3 );
    kitty::create_from_words( tt_maj, &_maj, &_maj + 1 );
    _e_storage->data.cache.insert( tt_maj );

    static uint64_t _ite = 0xd8;
    kitty::dynamic_truth_table tt_ite( 3 );
    kitty::create_from_words( tt_ite, &_ite, &_ite + 1 );
    _e_storage->data.cache.insert( tt_ite );

    static uint64_t _xor3 = 0x96;
    kitty::dynamic_truth_table tt_xor3( 3 );
    kitty::create_from_words( tt_xor3, &_xor3, &_xor3 + 1 );
    _e_storage->data.cache.insert( tt_xor3 );

    /* truth tables for constants */
    _e_storage->nodes[0].func = 0;

  }
  #pragma endregion constructors

#pragma region linking
  rig_network::e_signal_t rig_network::get_e_signal( rig_network::i_signal_t const& f )
  {
    return _e_storage->nodes[f].twin;
  }

  rig_network::i_signal_t rig_network::get_i_signal( rig_network::e_signal_t const& f )
  {
    return _i_storage->nodes[f].twin;
  }
#pragma endregion linking

#pragma region Primary I / O and constants
  rig_network::signal rig_network::get_constant( bool value = false )
  {
    return { 0, static_cast<uint64_t>( value ? 1 : 0 ) };
  }

  rig_network::signal rig_network::create_pi()
  {
    const auto e_index = _e_storage->get_index();
    auto& e_node = _e_storage->nodes.emplace_back();
    e_node.children.emplace_back( static_cast<uint64_t>(_e_storage->inputs.size()) );
    _e_storage->inputs.push_back( e_index );
    _e_storage->nodes[e_index].func = e_func_t::e_PI;

    const auto i_index = _i_storage->get_index();
    auto& i_node = _i_storage->nodes.emplace_back();
    i_node.children[0].data = i_node.children[1].data = _i_storage->inputs.size();
    _i_storage->inputs.emplace_back( i_index );
    _i_storage->nodes[i_index].func = i_func_t::i_PI;

    _e_storage->nodes[e_index].twin = { i_index, 0 };
    _i_storage->nodes[i_index].twin = { e_index, 0 };
    
    return { e_index, 0 };
  }

  uint32_t rig_network::create_po( rig_network::signal const& e_signal )
  {
    /* increase ref-count to children */
    _e_storage->nodes[e_signal.index].nfos++;
    auto const e_po_index = _e_storage->outputs.size();
    _e_storage->outputs.emplace_back( e_signal );
    
    i_signal_t i_signal = get_i_signal( e_signal );
    _i_storage->nodes[i_signal.index].nfos++;
    i_signal_t i_out_signal { i_signal.index, i_signal.complement ^ e_signal.complement };
    _i_storage->outputs.emplace_back( i_out_signal );

    return e_po_index;
  }

  bool rig_network::is_combinational() const
  {
    return true;
  }

  bool rig_network::is_constant( node const& n ) const
  {
    return n == 0;
  }

  bool rig_network::is_ci( node const& n ) const
  {
    return (_e_storage->nodes[n].children.size() == 1) && (_e_storage->nodes[n].children[0].index < num_pis() );
  }

  bool rig_network::is_pi( node const& n ) const
  {
    return (_e_storage->nodes[n].children.size() == 1) && (_e_storage->nodes[n].children[0].index < num_pis() );
  }
#pragma endregion Primary I / O and constants

#pragma region nodes and signals
  rig_network::node rig_network::get_node( rig_network::signal const& f ) const
  {
    return f.index;
  }

  rig_network::signal rig_network::make_signal( rig_network::node const& n ) const
  {
    return signal( n, 0 );
  }

  bool rig_network::is_complemented( rig_network::signal const& f ) const
  {
    return f.complement;
  }

  uint32_t rig_network::node_to_index( rig_network::node const& n ) const
  {
    return static_cast<uint32_t>( n );
  }

  rig_network::node rig_network::index_to_node( uint32_t index ) const
  {
    return index;
  }

  rig_network::node rig_network::ci_at( uint32_t index ) const
  {
    assert( index < _e_storage->inputs.size() );
    return *( _e_storage->inputs.begin() + index );
  }

  rig_network::signal rig_network::co_at( uint32_t index ) const
  {
    assert( index < _e_storage->outputs.size() );
    return *( _e_storage->outputs.begin() + index );
  }

  rig_network::node rig_network::pi_at( uint32_t index ) const
  {
    assert( index < _e_storage->inputs.size() );
    return *( _e_storage->inputs.begin() + index );
  }

  rig_network::signal rig_network::po_at( uint32_t index ) const
  {
    assert( index < _e_storage->outputs.size() );
    return *( _e_storage->outputs.begin() + index );
  }

  uint32_t rig_network::ci_index( rig_network::node const& n ) const
  {
    assert( _e_storage->nodes[n].children[0].data == _e_storage->nodes[n].children[1].data );
    return static_cast<uint32_t>( _e_storage->nodes[n].children[0].data );
  }

  uint32_t rig_network::pi_index( node const& n ) const
  {
    return static_cast<uint32_t>( _e_storage->nodes[n].children[0].data );
  }
#pragma endregion nodes and signals

#pragma region node and signal iterators
  template<typename Fn>
  void rig_network::foreach_node( Fn&& fn ) const
  {
    auto r = range<uint64_t>( _e_storage->nodes.size() );
    detail::foreach_element( r.begin(), r.end(), fn );
  }

  template<typename Fn>
  void rig_network::foreach_ci( Fn&& fn ) const
  {
    detail::foreach_element( _e_storage->inputs.begin(), _e_storage->inputs.end(), fn );
  }

  template<typename Fn>
  void rig_network::foreach_co( Fn&& fn ) const
  {
    using IteratorType = decltype( _e_storage->outputs.begin() );
    detail::foreach_element_transform<IteratorType, uint32_t>(
        _e_storage->outputs.begin(), _e_storage->outputs.end(), []( auto o ) { return o.data; }, fn );
  }

  template<typename Fn>
  void rig_network::foreach_pi( Fn&& fn ) const
  {
    detail::foreach_element( _e_storage->inputs.begin(), _e_storage->inputs.end(), fn );
  }

  template<typename Fn>
  void rig_network::foreach_po( Fn&& fn ) const
  {
    using IteratorType = decltype( _e_storage->outputs.begin() );
    detail::foreach_element_transform<IteratorType, uint32_t>(
        _e_storage->outputs.begin(), _e_storage->outputs.end(), []( auto o ) { return o.data; }, fn );
  }

  template<typename Fn>
  void rig_network::foreach_gate( Fn&& fn ) const
  {
    auto r = range<uint64_t>( 2u, _e_storage->nodes.size() ); /* start from 2 to avoid constants */
    detail::foreach_element_if(
        r.begin(), r.end(),
        [this]( auto n ) { return !is_ci( n ); },
        fn );
  }

  template<typename Fn>
  void rig_network::foreach_fanin( node const& n, Fn&& fn ) const
  {
    if ( n == 0 || is_ci( n ) )
      return;

    using IteratorType = decltype( _e_storage->outputs.begin() );
    detail::foreach_element_transform<IteratorType, uint32_t>(
        _e_storage->nodes[n].children.begin(), _e_storage->nodes[n].children.end(), []( auto f ) { return f.index; }, fn );
  }
#pragma endregion node and signal iterators

#pragma region unary functions
  rig_network::signal rig_network::create_buf( signal const& f )
  {
  //  auto [e_signal, is_new] = _create_known_node( { f }, e_func_t::e_BUF );
  //  _e_storage->nodes[e_signal.index].twin = _e_storage->nodes[f.index].twin;
  // return e_signal;
    return f;
  }

  rig_network::signal rig_network::create_not( signal const& f )
  {
    //auto [e_signal, is_new] = _create_known_node( { f }, e_func_t::e_BUF ^ 0x1 );
    //_e_storage->nodes[e_signal.index].twin = _e_storage->nodes[f.index].twin;
    //return e_signal;
    return !f;
  }

  bool rig_network::is_buf( node const& n )
  {
    _e_storage->nodes[n].func == e_func_t::e_BUF;
  }

  bool rig_network::is_not( node const& n )
  {
    _e_storage->nodes[n].func == (e_func_t::e_BUF ^ 0x1);
  }
#pragma endregion unary functions

rig_network::i_signal_t rig_network::i_create_and( i_signal_t a, i_signal_t b )
{
  /* order inputs */
  if ( a.index > b.index )
  {
    std::swap( a, b );
  }

  /* trivial cases */
  if ( a.index == b.index )
  {
    return ( a.complement == b.complement ) ? a : get_constant( false );
  }
  else if ( a.index == 0 )
  {
    return a.complement ? b : get_constant( false );
  }

  std::shared_ptr<i_storage_t>::element_type::node_type node;
  node.children[0] = a;
  node.children[1] = b;

  /* structural hashing */
  const auto it = _i_storage->hash.find( node );
  if ( it != _i_storage->hash.end() )
  {
    //assert( !is_dead( it->second ) );
    return { it->second, 0 };
  }

  const auto index = _i_storage->get_index();

  if ( index >= .9 * _i_storage->nodes.capacity() )
  {
    _i_storage->nodes.reserve( static_cast<uint64_t>( 3.1415f * index ) );
    _i_storage->hash.reserve( static_cast<uint64_t>( 3.1415f * index ) );
  }

  _i_storage->nodes.push_back( node );

  _i_storage->hash[node] = index;

  /* increase ref-count to children */
  _i_storage->nodes[a.index].nfos++;
  _i_storage->nodes[b.index].nfos++;

  return { index, 0 };
}

rig_network::i_signal_t rig_network::i_create_xor( i_signal_t a, i_signal_t b )
{
  i_signal_t f1 = i_create_and( a, !b );
  i_signal_t f2 = i_create_and( !a, b );
  return !i_create_and( !f1, !f2 );
}

#pragma region binary functions
  rig_network::signal rig_network::create_and( signal a, signal b )
  {
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }
    auto [e_signal, is_new] = _create_known_node( { a, b }, e_func_t::e_AND );
    if( is_new )
    {
      i_signal_t twin_a = _e_storage->nodes[get_node( a )].twin;
      i_signal_t twin_b = _e_storage->nodes[get_node( b )].twin;
      i_signal_t i_signal = i_create_and( { twin_a.index, twin_a.complement ^ a.complement }, { twin_b.index, twin_b.complement ^ b.complement } );
      _i_storage->nodes[i_signal.index].twin = e_signal;
      _e_storage->nodes[e_signal.index].twin = i_signal;
    }
    return e_signal;
  }

  rig_network::signal rig_network::create_nand( signal a, signal b )
  {
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }
    auto [e_signal, is_new] = _create_known_node( { a, b }, (e_func_t::e_AND ^ 1u ) );
    if( is_new )
    {
      i_signal_t twin_a = _e_storage->nodes[get_node( a )].twin;
      i_signal_t twin_b = _e_storage->nodes[get_node( b )].twin;
      i_signal_t i_signal = !i_create_and( { twin_a.index, twin_a.complement ^ a.complement }, { twin_b.index, twin_b.complement ^ b.complement } );
      _i_storage->nodes[i_signal.index].twin = e_signal;
      _e_storage->nodes[e_signal.index].twin = i_signal;
    }
    return e_signal;
  }

  rig_network::signal rig_network::create_or( signal a, signal b )
  {
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }
    auto [e_signal, is_new] = _create_known_node( { a, b }, e_func_t::e_OR );
    if( is_new )
    {
      i_signal_t twin_a = _e_storage->nodes[get_node( a )].twin;
      i_signal_t twin_b = _e_storage->nodes[get_node( b )].twin;
      i_signal_t i_signal = !i_create_and( { twin_a.index, 1u ^ twin_a.complement ^ a.complement }, { twin_b.index, 1u ^ twin_b.complement ^ b.complement } );
      _i_storage->nodes[i_signal.index].twin = e_signal;
      _e_storage->nodes[e_signal.index].twin = i_signal;
    }
    return e_signal;
  }

  rig_network::signal rig_network::create_nor( signal a, signal b )
  {
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }
    auto [e_signal, is_new] = _create_known_node( { a, b }, e_func_t::e_OR ^ 1u );
    if( is_new )
    {
      i_signal_t twin_a = _e_storage->nodes[get_node( a )].twin;
      i_signal_t twin_b = _e_storage->nodes[get_node( b )].twin;
      i_signal_t i_signal = i_create_and( { twin_a.index, 1u ^ twin_a.complement ^ a.complement }, { twin_b.index, 1u ^ twin_b.complement ^ b.complement } );
      _i_storage->nodes[i_signal.index].twin = e_signal;
      _e_storage->nodes[e_signal.index].twin = i_signal;
    }
    return e_signal;
  }

  rig_network::signal rig_network::create_lt( signal a, signal b )
  {
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }
    auto [e_signal, is_new] = _create_known_node( { a, b }, e_func_t::e_LT );
    if( is_new )
    {
      i_signal_t twin_a = _e_storage->nodes[get_node( a )].twin;
      i_signal_t twin_b = _e_storage->nodes[get_node( b )].twin;
      i_signal_t i_signal = i_create_and( { twin_a.index, 1u ^ twin_a.complement ^ a.complement }, { twin_b.index, twin_b.complement ^ b.complement } );
      _i_storage->nodes[i_signal.index].twin = e_signal;
      _e_storage->nodes[e_signal.index].twin = i_signal;
    }
    return e_signal;
  }

  rig_network::signal rig_network::create_ge( signal a, signal b )
  {
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }
    auto [e_signal, is_new] = _create_known_node( { a, b }, e_func_t::e_LT ^ 1u );
    if( is_new )
    {
      i_signal_t twin_a = _e_storage->nodes[get_node( a )].twin;
      i_signal_t twin_b = _e_storage->nodes[get_node( b )].twin;
      i_signal_t i_signal = !i_create_and( { twin_a.index, 1u ^ twin_a.complement ^ a.complement }, { twin_b.index, twin_b.complement ^ b.complement } );
      _i_storage->nodes[i_signal.index].twin = e_signal;
      _e_storage->nodes[e_signal.index].twin = i_signal;
    }
    return e_signal;
  }

  rig_network::signal rig_network::create_gt( signal a, signal b )
  {
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }
    auto [e_signal, is_new] = _create_known_node( { a, b }, e_func_t::e_GT );
    if( is_new )
    {
      i_signal_t twin_a = _e_storage->nodes[get_node( a )].twin;
      i_signal_t twin_b = _e_storage->nodes[get_node( b )].twin;
      i_signal_t i_signal = i_create_and( { twin_a.index, twin_a.complement ^ a.complement }, { twin_b.index, 1u ^ twin_b.complement ^ b.complement } );
      _i_storage->nodes[i_signal.index].twin = e_signal;
      _e_storage->nodes[e_signal.index].twin = i_signal;
    }
    return e_signal;
  }

  rig_network::signal rig_network::create_le( signal a, signal b )
  {
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }
    auto [e_signal, is_new] = _create_known_node( { a, b }, e_func_t::e_GT ^ 1u );
    if( is_new )
    {
      i_signal_t twin_a = _e_storage->nodes[get_node( a )].twin;
      i_signal_t twin_b = _e_storage->nodes[get_node( b )].twin;
      i_signal_t i_signal = !i_create_and( { twin_a.index, twin_a.complement ^ a.complement }, { twin_b.index, 1u ^ twin_b.complement ^ b.complement } );
      _i_storage->nodes[i_signal.index].twin = e_signal;
      _e_storage->nodes[e_signal.index].twin = i_signal;
    }
    return e_signal;
  }

  rig_network::signal rig_network::create_xor( signal a, signal b )
  {
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }
    auto [e_signal, is_new] = _create_known_node( { a, b }, e_func_t::e_XOR );
    i_signal_t twin_a = _e_storage->nodes[get_node( a )].twin;
    i_signal_t twin_b = _e_storage->nodes[get_node( b )].twin;
    i_signal_t i_signal = i_create_xor( { twin_a.index, twin_a.complement ^ a.complement }, { twin_b.index, twin_b.complement ^ b.complement } );
    _i_storage->nodes[i_signal.index].twin = e_signal;
    _e_storage->nodes[e_signal.index].twin = i_signal;
    return e_signal;
  }

  rig_network::signal rig_network::create_xnor( signal a, signal b )
  {
    if ( a.index > b.index )
    {
      std::swap( a, b );
    }
    auto [e_signal, is_new] = _create_known_node( { a, b }, e_func_t::e_XOR ^ 1u );
    i_signal_t twin_a = _e_storage->nodes[get_node( a )].twin;
    i_signal_t twin_b = _e_storage->nodes[get_node( b )].twin;
    i_signal_t i_signal = !i_create_xor( { twin_a.index, twin_a.complement ^ a.complement }, { twin_b.index, twin_b.complement ^ b.complement } );
    _i_storage->nodes[i_signal.index].twin = e_signal;
    _e_storage->nodes[e_signal.index].twin = i_signal;
    return e_signal;
  }

  bool rig_network::is_and( node const& n )
  {
    return _e_storage->nodes[n].func == e_AND;
  }
  bool rig_network::is_nand( node const& n )
  {
    return ( _e_storage->nodes[n].func == ( e_AND ^ 0x1 ) );
  }
  bool rig_network::is_or( node const& n )
  {
    return _e_storage->nodes[n].func == e_OR;
  }

  bool rig_network::is_nor( node const& n )
  {
    return _e_storage->nodes[n].func == ( e_OR ^ 0x1 );
  }

  bool rig_network::is_lt( node const& n )
  {
    return _e_storage->nodes[n].func == e_LT;
  }

  bool rig_network::is_ge( node const& n )
  {
    return _e_storage->nodes[n].func == ( e_LT ^ 0x1 );
  }

  bool rig_network::is_gt( node const& n )
  {
    return _e_storage->nodes[n].func == e_GT;
  }

  bool rig_network::is_le( node const& n )
  {
    return _e_storage->nodes[n].func == ( e_GT ^ 0x1 );
  }

  bool rig_network::is_xor( node const& n )
  {
    return _e_storage->nodes[n].func == e_XOR;
  }

  bool rig_network::is_xnor( node const& n )
  {
    return _e_storage->nodes[n].func == ( e_XOR ^ 0x1 );
  }
#pragma endregion binary functions

#pragma region arbitrary function
  std::tuple<rig_network::signal, bool> rig_network::_create_known_node( std::vector<signal> const& children, uint32_t literal )
  {
    std::shared_ptr<e_storage_t>::element_type::node_type node;
    std::copy( children.begin(), children.end(), std::back_inserter( node.children ) );
    node.func = literal;
    
    const auto it = _e_storage->hash.find( node );

    if ( it != _e_storage->hash.end() )
    {
      return std::make_pair( rig_network::signal{it->second, 0} , false );
    }

    const auto e_index = _e_storage->get_index();
    _e_storage->nodes.push_back( node );
    _e_storage->hash[node] = e_index;

    /* increase ref-count to children */
    for ( auto c : children )
    {
      _e_storage->nodes[c.index].nfos++;
    }

    return std::make_pair<signal, bool>( rig_network::signal{ e_index, 0 }, true );
  }

  rig_network::signal rig_network::_create_node( std::vector<signal> const& children, uint32_t literal )
  {
    std::shared_ptr<e_storage_t>::element_type::node_type node;
    std::copy( children.begin(), children.end(), std::back_inserter( node.children ) );
    node.func = literal;

    const auto it = _e_storage->hash.find( node );
    if ( it != _e_storage->hash.end() )
    {
      return it->second;
    }

    const auto e_index = _e_storage->get_index();
    _e_storage->nodes.push_back( node );
    _e_storage->hash[node] = e_index;

    /* increase ref-count to children */
    for ( auto c : children )
    {
      _e_storage->nodes[c].nfos++;
    }

    // synthesize
    std::array<i_signal_t, 4u> i_children {};
    for( auto i{0}; children.size(); ++i )
      i_children[i] = _e_storage->nodes[children[i].index].twin;
    auto i_signal = synthesize_twin( i_children, literal );
    _e_storage->nodes[e_index].twin = i_signal;
    _i_storage->nodes[i_signal.index].twin = { e_index, 0 };

    return { e_index, 0 };
  }

  bool rig_network::is_function( node const& n ) const
  {
    return n > 0 && !is_ci( n );
  }

  rig_network::i_signal_t rig_network::synthesize_twin( std::array<i_signal_t, 4u> & i_children, uint32_t literal )
  {
    kitty::dynamic_truth_table const& tt = _e_storage->data.cache[literal];
    return synthesize_twin_rec( i_children, tt );
  }

  rig_network::i_signal_t rig_network::synthesize_twin_rec( std::array<i_signal_t, 4u> & i_children, kitty::dynamic_truth_table const& tt )
  {
    if( i_children.size() < 4u )
    {
      return boolean_matching( i_children, tt );
    }
    return synthesize_twin_rec( i_children, tt );
  }

  rig_network::i_signal_t rig_network::boolean_matching( std::array<i_signal_t, 4u> & i_children, kitty::dynamic_truth_table const& tt )
  {
    return {0,0};
//    kitty::static_truth_table<4u> tt_s = kitty::extend_to<4u>( tt );
//
//    auto [func_npn, neg, perm] = exact_npn_canonization( tt_s );
//    auto const structures = _database.get_supergates( func_npn );
//    bool phase = ( neg >> 4 == 1 ) ? true : false;
//
//    for( auto i{0}; i<i_children.size(); ++i )
//    {
//      if( ( neg >> i ) & 0x1 == 0x1 )
//        i_children[i] = !i_children[i];
//    }
//    std::array<i_signal_t, 4> leaves;
//    for( auto i{0}; i<4; ++i )
//    {
//      leaves[i] = i_children[perm[i]];
//    }
//
//    auto & db = _database.get_database();
//    i_signal_t i_signal = {0,0};//create_twin_network( db.get_node( structures->at(0).root ), leaves );
//    
//    return phase != db.is_complemented( structures->at(0).root ) ? !i_signal : i_signal;

  }
#pragma endregion arbitrary function

#pragma region structural properties
  size_t rig_network::size() const
  {
    return static_cast<uint32_t>( _e_storage->nodes.size() );
  }

  size_t rig_network::num_cis() const
  {
    return static_cast<uint32_t>( _e_storage->inputs.size() );
  }

  size_t rig_network::num_cos() const
  {
    return static_cast<uint32_t>( _e_storage->outputs.size() );
  }

  size_t rig_network::num_pis() const
  {
    return static_cast<uint32_t>( _e_storage->inputs.size() );
  }

  size_t rig_network::num_pos() const
  {
    return static_cast<uint32_t>( _e_storage->outputs.size() );
  }

  size_t rig_network::num_gates() const
  {
    return static_cast<uint32_t>( _e_storage->nodes.size() - _e_storage->inputs.size() - 1 );
  }

  size_t rig_network::fanin_size( node const& n ) const
  {
    return static_cast<uint32_t>( _e_storage->nodes[n].children.size() );
  }

  size_t rig_network::fanout_size( node const& n ) const
  {
    return _e_storage->nodes[n].nfos & UINT32_C( 0x7FFFFFFF );
  }

  size_t rig_network::incr_fanout_size( node const& n ) const
  {
    return _e_storage->nodes[n].nfos++ & UINT32_C( 0x7FFFFFFF );
  }

  size_t rig_network::decr_fanout_size( node const& n ) const
  {
    return --_e_storage->nodes[n].nfos & UINT32_C( 0x7FFFFFFF );
  }

#pragma endregion structural properties

#pragma region functional properties
  kitty::dynamic_truth_table rig_network::node_function( const node& n ) const
  {
    return _e_storage->data.cache[_e_storage->nodes[n].func];
  }
#pragma endregion functional properties

} // namespace rils

} // namespace mockturtle
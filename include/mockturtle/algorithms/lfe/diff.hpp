diff --git a/include/mockturtle/algorithms/lfe/sim_decomposition_fastS.hpp b/include/mockturtle/algorithms/lfe/sim_decomposition_fastS.hpp
index 1190bd73..62dd5b62 100644
--- a/include/mockturtle/algorithms/lfe/sim_decomposition_fastS.hpp
+++ b/include/mockturtle/algorithms/lfe/sim_decomposition_fastS.hpp
@@ -60,7 +60,7 @@ namespace mockturtle
 /*! \brief Parameters for sim_decomposition_fastS algorithm */
 struct sim_decomposition_fastS_params
 {
-  bool verbose{true};
+  bool verbose{false};
   uint32_t max_sup{2};
   bool is_informed{true};
   bool is_size_aware{false};
@@ -166,7 +166,6 @@ namespace detail
           kitty::print_binary(chj_res.dtt);
           std::cout << std::endl;
         }
-
         return fc;
       }
 
@@ -178,9 +177,6 @@ namespace detail
         bool is_success = false;
         TT on_xi;
         TT off_xi;
-        double Llim = 1.00;//support.size() < 256 ? 0.95 : 1.00;
-        double Rlim = 1.00;//support.size() < 256 ? 1.05 : 1.00;
-        double RTIO = 0.00;//support.size() < 256 ? 0.01 : 0.0;
 
         chj_result chj_new_node;
         /*{
@@ -228,17 +224,10 @@ namespace detail
           for( uint32_t j = i+1; j < sorted_indeces.size(); ++j )
           {
             support_pat_pointers[1] = &X[support[sorted_indeces[j]]].pat;
-            //if( vect_I[j] < vect_I[i] )
-
-            double rtio = vect_I[j] == 0 ? 0 : abs( (vect_I[i] - vect_I[j] )/vect_I[i]);
-            //if( vect_I[j] >= vect_I[i] )
-            //  std::cout << rtio << std::endl;
-            if( ( (vect_I[i] > vect_I[j]) && ( rtio > RTIO ) )  )
+            if( vect_I[i] != vect_I[j] )
               break;
             else
             {
-              //\std::cout << i << " " << vect_I[i] << " " << j << " " << vect_I[j] << std::endl; 
-
               create_candidates_result<TT> candidates = create_candidates_method( support_pat_pointers, &on_f );
               for( uint32_t k = 0; k < candidates.dtt_v.size(); ++k )
               {
@@ -246,35 +235,24 @@ namespace detail
                 off_xi = amask & ~candidates.pat_v[k];
                 Inew = information( on_xi, off_xi, on_f, off_f );
 ////////////////////////////////////////////////////////////////
+                TT xl = amask & X[support[sorted_indeces[i]]].pat;
+                TT xr = amask & X[support[sorted_indeces[j]]].pat;
+                TT ym = amask & on_f;
+                TT xn = on_xi;
 
+                std::vector<TT*> Xptr;
+                Xptr = std::vector{&xl, &xr};
+                double Iij = kitty::mutual_information( Xptr, &ym );
+                Xptr = std::vector{&xn};
+                double In = kitty::mutual_information( Xptr, &ym );
+                Xptr = std::vector{&xl,&xn};
+                double Iin = kitty::mutual_information( Xptr, &ym );
+                Xptr = std::vector{&xr, &xn};
+                double Ijn = kitty::mutual_information( Xptr, &ym );
+                Xptr = std::vector{&xl, &xr, &xn};
+                double Iijn = kitty::mutual_information( Xptr, &ym );
 
-                if( Inew > Imax )
-                {
-                  TT xl = amask & X[support[sorted_indeces[i]]].pat;
-                  TT xr = amask & X[support[sorted_indeces[j]]].pat;
-                  TT ym = amask & on_f;
-                  TT xn = on_xi;
-
-                  std::vector<TT*> Xptr;
-                  //A Xptr = std::vector{&xl, &xn};
-                  //A double Iin = kitty::mutual_information( Xptr, &ym );
-                  //A Xptr = std::vector{&xr, &xn};
-                  //A double Ijn = kitty::mutual_information( Xptr, &ym );
-                  //A if( Iin == Ijn )
-                  {
-                  Xptr = std::vector{&xl, &xr};
-                  double Iij = kitty::mutual_information( Xptr, &ym );
-                  //A if( Iij == Iin )
-                  {
-                  Xptr = std::vector{&xn};
-                  double In = kitty::mutual_information( Xptr, &ym );
-                  Xptr = std::vector{&xl,&xn};
-                  if( ( In >= Llim*Iij && In <= Rlim*Iij  ) )
-                  {
-                  //  Xptr = std::vector{&xl, &xr, &xn};
-                  //  double Iijn = kitty::mutual_information( Xptr, &ym );
-                 //if( Iijn == In && Iin == In && Ijn == In && Iij == In )
-                 //if(  Iijn == Iin )//&& Iin == In && Ijn == In && Iij == In )
+                if( Inew > Imax && Iijn == In && Iin == In && Ijn == In && Iij == In )
 ///////////////////////////////////////////////////////
                 {
                   std::vector<signal<Ntk>> children;
@@ -294,10 +272,6 @@ namespace detail
                     is_success = true;
                   }
                 }
-                }
-                  }
-                  }
-                }
 ///////////////////////////////////////////////////////
               }
             }
@@ -313,8 +287,8 @@ namespace detail
           //std::cout << support.size() << std::endl;
           X.push_back( ntk.sim_patterns[ ntk.get_node_pattern( fc ) ] );
           std::cout << children[0] << " " << children[1] << " " << chj_new_node.tt << std::endl;  
-          //support.erase( support.begin() + std::max<uint32_t>( new_node_support_indeces[0],new_node_support_indeces[1] ) );
-          //support.erase( support.begin() + std::min<uint32_t>( new_node_support_indeces[0],new_node_support_indeces[1] ) );
+          support.erase( support.begin() + std::max<uint32_t>( new_node_support_indeces[0],new_node_support_indeces[1] ) );
+          support.erase( support.begin() + std::min<uint32_t>( new_node_support_indeces[0],new_node_support_indeces[1] ) );
         }
         return is_success;
       }
@@ -445,26 +419,13 @@ namespace detail
     
       signal<Ntk> idsd_step( std::vector<uint32_t> support, TT amask, TT xmask, bool branch_on_last = false )
       {
-        uint32_t BFcnt=0;
-        std::cout << "b1 " << std::endl; 
-
-        //for( auto x : support )
-        //  std:: cout << x << " "; 
-        std::cout << support.size() << std::endl;
-        std::cout << "d" << BFcnt++ << std::endl;
-
-        //std::cout << support.size() << std::endl;
-        //std::cout << num_nobranch << std::endl;
         uint32_t n_ones = kitty::count_ones(amask);
 
         if( n_ones == 0 )
-        {
-          return ntk.get_constant( false ); 
-        }   
+          return ntk.get_constant( false );    
         if( support.size() == 0 )
-        {
           return ntk.get_constant( false );
-        }
+
         TT on_f(n_bits); 
         TT off_f(n_bits); 
 
@@ -472,20 +433,15 @@ namespace detail
         off_f = amask&(~(xmask^Y.pat ));
 
         if( kitty::count_ones( on_f ) == 0 ) // contraddiction
-        {
           return ntk.get_constant( false );
-        }
         else if( kitty::count_ones( on_f ) == n_ones ) // tautology
-        {
           return ntk.get_constant( true );
-        }
 
         uint32_t bidx = 0;
         signal<Ntk> bsig;
         double Inew = -std::numeric_limits<double>::max();
         double Imax = -std::numeric_limits<double>::max();
         uint32_t max_fanin_size = std::numeric_limits<uint32_t>::max();
-        
 
         std::vector<uint32_t> to_be_deleted;
         TT on_x(n_bits);
@@ -503,8 +459,6 @@ namespace detail
             off_x = amask & ~X[support[bidx]].pat;
             bsig = X[support[bidx]].sig;
 
-            Imax = information( on_x, off_x, on_f, off_f );
-
             if( on_x == on_f )
               return X[support[bidx]].sig;
             if( on_x == off_f )
@@ -559,9 +513,7 @@ namespace detail
             off_xi = amask & ~X[support[i]].pat;
 
             if( on_xi == on_f )
-            {
               return X[support[i]].sig;
-            }
             if( on_xi == off_f )
             {
               signal<klut_network>fo = ntk.create_not( X[support[i]].sig );
@@ -576,24 +528,18 @@ namespace detail
           bidx = 0;
         }
 
-
+          
         if( to_be_deleted.size() > 0 )
         {
-          std::cout << "H" << BFcnt++ << std::endl;
           std::reverse(to_be_deleted.begin(), to_be_deleted.end());
           uint32_t count = 0;
           for( uint32_t i = 0; i < to_be_deleted.size() ; ++i )
           {
-            std::cout << support.size() << " " << to_be_deleted[i] << std::endl;
-            assert( support.size() > to_be_deleted[i] );
-            std::cout << "delete " << to_be_deleted[i] << " from " << support.size() << std::endl;
             support.erase( support.begin()+to_be_deleted[i] ) ;
             if( to_be_deleted[i] < bidx )
               count++;
           }
-          bidx = uint32_t(bidx-count);
-
-          std::cout << "bidx " << bidx << std::endl;
+          bidx -= count;
         }
 
         if( support.size() == 0 )
@@ -621,7 +567,6 @@ namespace detail
         std::vector<uint32_t> reduced_support = support;
 
         reduced_support.erase( reduced_support.begin() + bidx );
-        std::cout << "A4" << std::endl;
 
         std::vector<uint32_t> pis_support;
         if( ps.try_xor )
@@ -692,15 +637,12 @@ namespace detail
           }
         }
 
-        if( !branch_on_last && ps.try_bottom_decomposition ) //support.size() < 256 && 
+        if( !branch_on_last && ps.try_bottom_decomposition )
         {
           if( ps.is_informed )
           {
             if ( try_bottom_decomposition( support, amask, on_f, off_f, Imax ) )
-            {
-              //num_nobranch += 1;
               return idsd_step( support, amask, xmask, true );
-            }
           }
           else
           { 
@@ -708,38 +650,19 @@ namespace detail
           }
         }
         if( ps.is_size_aware ) 
-        {
           clear_fanin_size( bsig );
-        }
-        std::cout << "A10-1" << std::endl;
-        std::cout << kitty::count_ones( amask0 ) << std::endl;
-        std::cout << kitty::count_ones( ~amask0 ) << std::endl;
-        std::cout << reduced_support.size() << std::endl;
-        //for( auto x : reduced_support )
-        //  std::cout << x << " ";
-        //std::cout << std::endl;
-        
-        klut_network::signal F0 = idsd_step( reduced_support, amask0, xmask0 );
-        std::cout << "A10+1" << std::endl;
 
-        klut_network::signal f0 = ntk.create_and( ntk.create_not( bsig ), F0 );
-        std::cout << "A11" << std::endl;
+        signal<klut_network> F0 = idsd_step( reduced_support, amask0, xmask0 );
+        signal<klut_network> f0 = ntk.create_and( ntk.create_not( bsig ), F0 );
 
         signal<klut_network> F1 = idsd_step( reduced_support, amask1, xmask1 );
-        std::cout << "A12" << std::endl;
-
         signal<klut_network> f1 = ntk.create_and( bsig, F1 );
-        std::cout << "A13" << std::endl;
-
 
         signal<klut_network> Fnew = ntk.create_or( f1, f0 );
 
         if( ps.verbose )
           std::cout << Fnew << "= ite(" << bsig << "," << F1 << "," << F0 << ")" << std::endl;
 
-        std::cout << "b0" << std::endl;
-
-
         return Fnew;
       }
 
@@ -749,17 +672,11 @@ namespace detail
 
         for( uint32_t i = 0; i < X.size(); ++i )
           support.push_back( i );
-        
-        if( X.size() > 256 )
-          size_thresh = uint32_t( 0.5*X.size() );
       
 
         TT xmask( Y.pat.num_bits() );
         TT amask = ~xmask;
-        std::cout << "a1" << std::endl;
         return idsd_step( support, amask, xmask );
-        std::cout << "a0" << std::endl;
-
       }
     
     private:
@@ -768,7 +685,6 @@ namespace detail
       kitty::partial_truth_table target;
       sim_pattern<Ntk> Y;
       uint32_t n_bits;
-      uint32_t size_thresh;
       bool branch_on_all;
     public:
       std::vector<double> Iactive;

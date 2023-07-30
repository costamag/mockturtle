#include <kitty/constructors.hpp>
#include <kitty/partial_truth_table.hpp>
#include <iostream>
#include <fstream>
#include<sstream>

std::vector<int> getNumberFromString(std::string s) {
    std::vector<int> res; 
   std::stringstream str_strm;
   str_strm << s; //convert the string s into stringstream
   std::string temp_str;
   int temp_int;
   while(!str_strm.eof()) {
      str_strm >> temp_str; //take words into temp_str one by one
      if(std::stringstream(temp_str) >> temp_int) { //try to convert string to int
        res.push_back( temp_int );
        //std::cout << temp_int << " ";
      }
      //temp_str = ""; //clear temp string
   }
   return res;
}

std::vector<int> greedy_set_covering( std::vector<kitty::partial_truth_table> S, std::vector<int> W, kitty::partial_truth_table T )
{
    assert( S.size() == W.size() );
    std::vector<int> res;
    std::vector<int> support;
    for( int i{0}; i<S.size(); ++i ) support.push_back( i );
    double prize;
    int candidate;
    while( kitty::count_ones(T) > 0 )
    {
        double best_prize = 10000000;
        for( int i{0}; i<support.size(); ++i )
        {
            int nNewlyAdded = kitty::count_ones( T & S[support[i]] );
            if( nNewlyAdded <= 0 ) continue;
            prize = W[support[i]]*1.0/( nNewlyAdded );
            //printf("%f\n", prize);
            if( prize < best_prize )
            {
                best_prize = prize;
                candidate = i;
            }
        }
        T &= ~S[support[candidate]];
        printf("chosen %d to go %d\n", candidate, kitty::count_ones(T));

        res.push_back(support[candidate]);
        support.erase( support.begin() + candidate );
    }
    return res;
}

std::vector<int> mod_greedy_set_covering( std::vector<kitty::partial_truth_table> S, std::vector<int> W, kitty::partial_truth_table T )
{
    assert( S.size() == W.size() );
    std::vector<int> res;
    std::vector<int> support;
    for( int i{0}; i<S.size(); ++i ) support.push_back( i );
    double prize;

    while( kitty::count_ones(T) > 0 )
    {
        kitty::partial_truth_table G0, G1, G2, G0t, G1t, G2t;
        double best_prize = 10000000;
        int n = support.size();
        int candidate=0;

        if( n > 2 )
        {
            G0 = S[support[n-2]];
            G1 = S[support[n-1]];
            G2 = S[support[n]];

            for( int i{0}; i < n; ++i )
            {
                int nNewlyAdded = kitty::count_ones( T & S[support[i]] );
                if( nNewlyAdded <= 0 ) continue;
                prize = W[support[i]]*1.0/( nNewlyAdded );
                //printf("%f\n", prize);
                if( prize < best_prize )
                {
                    best_prize = prize;
                    candidate = i;
                }
                if( n-i >= 3 )
                {
                    int nn = n-i-3;
                    G0t = S[support[nn-2]];
                    G1t = G1 | G2;
                    G2t = G0 | ( G1 & G2 ); 
                    G0 = G0t;
                    G1 = G1t;
                    G2 = G2t;
                }
            }
            G0 = T&(G1 ^ G2);

        }
        else if( n == 2 )
        {
            G1 = S[support[n-2]];
            G2 = S[support[n-1]];   
            G0 = T&(G1 ^ G2);
        }
        else
        {
            G0 = T&S[support[0]];
        }

        if( kitty::count_ones(G0) == 0 || n == 1 )
        {
            T &= ~S[support[candidate]];
            printf("chosen %d to go %d : essentials = %d\n", candidate, kitty::count_ones(T), kitty::count_ones(G0));
            res.push_back(support[candidate]);
            support.erase( support.begin() + candidate );
        }
        else
        {
            std::vector<int> to_erase;
            for( int i{0}; i < n; ++i )
            {
                if( kitty::count_ones( G0 & S[support[i]] & T ) > 0 )
                {
                    T &= ~S[support[i]];
                    printf("chosen %d to go %d : essentials = %d\n", support[i], kitty::count_ones(T), kitty::count_ones(G0));
                    res.push_back(support[i]);
                    to_erase.push_back(i);
                }
            }
            for( int iEr{int(to_erase.size()-1)}; iEr >= 0; --iEr )
                support.erase( support.begin() + to_erase[iEr] );
        }
    }
    return res;
}

std::vector<int> find_essential( std::vector<kitty::partial_truth_table> M, kitty::partial_truth_table T, std::vector<int> support )
{
    std::vector<kitty::partial_truth_table> Gs;
    for( int i{0}; i<support.size(); ++i )
        Gs.push_back( M[support[i]] & T );
    kitty::partial_truth_table G1, G2; 
    if( support.size() > 1 )
    {
        for( int n{int(support.size()-1)}; n >= 2; --n )
        {
            G1 = Gs[n] | Gs[n-1];
            G2 = Gs[n-2] | ( Gs[n] & Gs[n-1] );
            Gs[n-1] = G1;
            Gs[n-2] = G2;
        }
        G1 = Gs[0] ^ Gs[1];
    }
    printf("nEssentials = %d\n", kitty::count_ones(G1));
    return support;   
}

int main( int argc, char * argv[] )
{
  std::vector<kitty::partial_truth_table> M;
  std::vector<kitty::partial_truth_table> Mt;
  int nRows;
  int nCols;

  assert( argc >= 2 );
  std::string file_name = "../experiments/set_covering/";
  file_name += argv[1];
  std::string line;
  std::ifstream myfile (file_name);
  bool is_first = true;
  int id{-1};
  std::vector<int> weights;
  std::vector<int> support;

  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      if( is_first )
      {
        //std::cout << line << '\n';
        is_first = false;
        auto specs = getNumberFromString(line);
        nRows = specs[0];
        nCols = specs[1];
        //std::cout << nRows << " " << nCols << std::endl;
        for( int iCol{0}; iCol < nCols; ++iCol )
            M.emplace_back( nRows );
        for( int iRow{0}; iRow < nRows; ++iRow )
            Mt.emplace_back( nCols );
      }
      else
      {
        auto specs = getNumberFromString(line);
        int w = specs[0];
        weights.push_back( w );
        support.push_back( id );
        
        int nR = specs[1];
        specs.erase( specs.begin() );
        specs.erase( specs.begin() );
        for( auto iRow : specs )
        {
            kitty::set_bit( M[id], iRow-1 );
            kitty::set_bit( Mt[iRow-1], id );
        }
      }
      id++;
    }
    myfile.close();
  }
  else std::cout << "Unable to open file " << file_name << std::endl;

kitty::partial_truth_table target(nRows);
    for( auto tt : M )
        target |= tt;
    printf("%d/%d", kitty::count_ones(target), target.num_bits());


    find_essential( M, target, support );
    auto RES = greedy_set_covering( M, weights, target );
    for( auto r : RES )
        printf("%d ", r );
    printf(" : |S|=%d\n", RES.size());

    auto RES2 = mod_greedy_set_covering( M, weights, target );
    for( auto r : RES2 )
        printf("%d ", r );
    printf(" : |S|=%d\n", RES2.size());

    return 0;
}
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
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file mnist_reader.hpp
  \brief handler of the mnist data

  \author Andrea Costamagna
*/
#pragma once

#include <kitty/partial_truth_table.hpp>
#include <kitty/bit_operations.hpp>
#include <kitty/operations.hpp>  
#include <kitty/operators.hpp>  
#include "ml_rng2.hpp"
#include <stdio.h>
#include <stack>
#include <iostream>
#include <string>
#include <fstream>
#include <bit>
#include <bitset>
#include <cstdint>
#include <fmt/format.h>
#include <random>
#include <set>


namespace mockturtle
{

namespace mcts
{

    using PTT = kitty::partial_truth_table;

    int reverseInt( int );
    
    int read_pixel( char );


    int reverseInt (int i) 
    {
        unsigned char c1, c2, c3, c4;

        c1 = i & 255;
        c2 = (i >> 8) & 255;
        c3 = (i >> 16) & 255;
        c4 = (i >> 24) & 255;

        return ((int)c1 << 24) + ((int)c2 << 16) + ((int)c3 << 8) + c4;
    }

    template<class IETYPE> int read_pixel( char );

    std::vector<PTT> read_mnist_image_bin( std::string fileName, int thr = 60000 )
    {
        std::vector<PTT> X;
        std::ifstream file ( fileName.c_str() );
        if (file.is_open())
        {
            int magic_number=0;
            int number_of_images=0;
            int n_rows=0;
            int n_cols=0;
            file.read((char*)&magic_number,sizeof(magic_number)); 
            magic_number= reverseInt(magic_number);
            file.read((char*)&number_of_images,sizeof(number_of_images));
            number_of_images= reverseInt(number_of_images);
            file.read((char*)&n_rows,sizeof(n_rows));
            n_rows= reverseInt(n_rows);
            file.read((char*)&n_cols,sizeof(n_cols));
            n_cols= reverseInt(n_cols);
            for(int i=0;i<n_rows*n_cols;++i)
                X.emplace_back();

            for(int i=0;i<number_of_images;++i)
            {
                if( i > thr ) break;
                for(int r=0;r<n_rows;++r)
                {
                    for(int c=0;c<n_cols;++c)
                    {
                        unsigned char temp=0;
                        file.read((char*)&temp,sizeof(temp));
                        X[r*28+c].add_bit( (int)temp > 1 );
                    }
                }
            }
        }
        else
        {
            printf("FILE NOT FOUND\n");
            assert(0);
        }
        return X;
    }

    std::vector<PTT> read_mnist_label_04_59( std::string fileName, int thr = 60000 )
    {
        std::vector<PTT> Y;
        std::ifstream file ( fileName.c_str() );
        if (file.is_open())
        {
            int magic_number=0;
            int number_of_images=0;
            int n_rows=0;
            int n_cols=0;
            file.read((char*)&magic_number,sizeof(magic_number)); 
            magic_number= reverseInt(magic_number);
            file.read((char*)&number_of_images,sizeof(number_of_images));
            number_of_images= reverseInt(number_of_images);
            Y.emplace_back();
            for(int i=0;i<number_of_images;++i)
            {
                if( i > thr ) break;
                unsigned char temp=0;
                file.read((char*)&temp,sizeof(temp));
                Y[0].add_bit( (int)(temp) > 4 );
                //if( i < 10 ) printf("%d %d\n", temp, (int)(temp) > 4 );
            }
        }
        else
        {
            printf("FILE NOT FOUND\n");
            assert(0);
        }
        return Y;
    }

    std::vector<int> read_mnist_label_same( std::string fileName, int thr = 60000 )
    {
        std::vector<int> Y;
        std::ifstream file ( fileName.c_str() );
        if (file.is_open())
        {
            int magic_number=0;
            int number_of_images=0;
            int n_rows=0;
            int n_cols=0;
            file.read((char*)&magic_number,sizeof(magic_number)); 
            magic_number= reverseInt(magic_number);
            file.read((char*)&number_of_images,sizeof(number_of_images));
            number_of_images= reverseInt(number_of_images);
            for(int i=0;i<number_of_images;++i)
            {
                if( i > thr ) break;
                unsigned char temp=0;
                file.read((char*)&temp,sizeof(temp));
                Y.push_back( (int)(temp) );
                //if( i < 10 ) printf("%d %d\n", temp, (int)(temp) > 4 );
            }
        }
        else
        {
            printf("FILE NOT FOUND\n");
            assert(0);
        }
        return Y;
    }


    std::vector<PTT> read_mnist_label_10( std::string fileName, int thr = 60000 )
    {
        std::vector<PTT> Y;
        std::ifstream file ( fileName.c_str() );
        if (file.is_open())
        {
            int magic_number=0;
            int number_of_images=0;
            int n_rows=0;
            int n_cols=0;
            file.read((char*)&magic_number,sizeof(magic_number)); 
            magic_number= reverseInt(magic_number);
            file.read((char*)&number_of_images,sizeof(number_of_images));
            number_of_images= reverseInt(number_of_images);
            for( int i{0}; i<5; ++i )
                Y.emplace_back();
            for(int i=0;i<number_of_images;++i)
            {
                if( i > thr ) break;
                unsigned char temp=0;
                file.read((char*)&temp,sizeof(temp));
                switch (temp)
                {
                    case 0: Y[0].add_bit(1); Y[1].add_bit(0); Y[2].add_bit(0); Y[3].add_bit(0); Y[4].add_bit(0);    break;
                    case 1: Y[0].add_bit(1); Y[1].add_bit(0); Y[2].add_bit(0); Y[3].add_bit(0); Y[4].add_bit(1);    break;
                    case 2: Y[0].add_bit(0); Y[1].add_bit(0); Y[2].add_bit(0); Y[3].add_bit(1); Y[4].add_bit(0);    break;
                    case 3: Y[0].add_bit(0); Y[1].add_bit(0); Y[2].add_bit(0); Y[3].add_bit(1); Y[4].add_bit(1);    break;
                    case 4: Y[0].add_bit(0); Y[1].add_bit(0); Y[2].add_bit(1); Y[3].add_bit(0); Y[4].add_bit(0);    break;
                    case 5: Y[0].add_bit(0); Y[1].add_bit(0); Y[2].add_bit(1); Y[3].add_bit(0); Y[4].add_bit(1);    break;
                    case 6: Y[0].add_bit(0); Y[1].add_bit(0); Y[2].add_bit(1); Y[3].add_bit(1); Y[4].add_bit(0);    break;
                    case 7: Y[0].add_bit(0); Y[1].add_bit(0); Y[2].add_bit(1); Y[3].add_bit(1); Y[4].add_bit(1);    break;
                    case 8: Y[0].add_bit(1); Y[1].add_bit(1); Y[2].add_bit(0); Y[3].add_bit(0); Y[4].add_bit(0);    break;
                    case 9: Y[0].add_bit(1); Y[1].add_bit(1); Y[2].add_bit(0); Y[3].add_bit(0); Y[4].add_bit(1);    break;
                default:
                    break;
                }
            }
        }
        else
        {
            printf("FILE NOT FOUND\n");
            assert(0);
        }
        return Y;
    }

    void print_mnist_image( std::vector<PTT> * pX, std::vector<PTT> * pY, int idx )
    {
        printf("IMAGE: ");
        for( int i{0}; i<pY->size(); ++i )
            printf("%d\n", kitty::get_bit( (*pY)[i], idx ));
        printf("\n");


        for( int r{0}; r<28; ++r )
        {
            for( int c{0}; c<28; ++c)
                printf("%d", kitty::get_bit( (*pX)[r*28+c], idx ));
            printf("\n");
        }
        printf("\n");
    }


    struct binary_classification_dataset_t
    {
        std::vector<PTT> x_train;
        std::vector<PTT> x_valid;
        std::vector<PTT> x_test;
        std::vector<std::vector<uint32_t>> v_x_train;
        std::vector<std::vector<uint32_t>> v_x_valid;
        std::vector<std::vector<uint32_t>> v_x_test;
        std::vector<uint32_t> v_y_train;
        std::vector<uint32_t> v_y_valid;
        std::vector<uint32_t> v_y_test;
        PTT y_train;
        PTT m_train;
        PTT y_valid;
        PTT m_valid;
        PTT y_test;
        PTT m_test;
    };

    binary_classification_dataset_t linearly_separable_dataset( int nBits, int nTrain, int nValid, int nTest )
    {
        assert(nBits < 32);
        binary_classification_dataset_t data;
        std::uniform_int_distribution<> distrib(0, ( 1u << nBits ) - 1 );

        PTT train0(nTrain);
        PTT valid0(nValid);
        PTT test0(nTest); 

        for( int i{0}; i<2*nBits; ++i )
        {
            data.x_train.push_back( train0 );
            data.x_valid.push_back( valid0 );
            data.x_test.push_back( test0 );
        }

        data.y_train = train0;
        data.y_valid = valid0;
        data.y_test = test0;
        data.m_train = ~train0;
        data.m_valid = ~valid0;
        data.m_test  = ~test0;

        int i=0;
        std::set<std::pair<uint32_t, uint32_t>> set;
        int rep=0;
        while( i<nTrain && rep < 1000 )
        {
            rep++;
            uint32_t x1 = distrib(ml_gen2);
            uint32_t x2 = distrib(ml_gen2);
            std::pair<uint32_t, uint32_t> pair = std::make_pair( x1, x2 );
            if( set.find( pair ) != set.end() )
                continue;
            set.insert(pair);
            i++;
            data.v_x_train.push_back( {x1, x2} );
            for( int j{0}; j<nBits; ++j )
            {
                if( ( x1 >> j ) & 1u == 1u )
                    kitty::set_bit( data.x_train[j], i );
                else
                    kitty::clear_bit( data.x_train[j], i );

                if( ( x2 >> j ) & 1u == 1u )
                    kitty::set_bit( data.x_train[j+nBits], i );
                else
                    kitty::clear_bit( data.x_train[j+nBits], i );
            }
            if( x2 > x1 )
                kitty::set_bit( data.y_train, i );
            else
                kitty::clear_bit( data.y_train, i );
            data.v_y_train.push_back( kitty::get_bit( data.y_train, i ) );
        }

        for( int i{0}; i<nValid; ++i )
        {
            uint32_t x1 = distrib(ml_gen2);
            uint32_t x2 = distrib(ml_gen2);
            data.v_x_valid.push_back( {x1, x2} );
            for( int j{0}; j<nBits; ++j )
            {
                if( ( x1 >> j ) & 1u == 1u )
                    kitty::set_bit( data.x_valid[j], i );
                else
                    kitty::clear_bit( data.x_valid[j], i );
                if( ( x2 >> j ) & 1u == 1u )
                    kitty::set_bit( data.x_valid[j+nBits], i );
            }
            if( x2 > x1 )
                kitty::set_bit( data.y_valid, i );
            else
                kitty::clear_bit( data.y_valid, i );
            data.v_y_valid.push_back( kitty::get_bit( data.y_valid, i ) );
        }

        for( int i{0}; i<nTest; ++i )
        {
            uint32_t x1 = distrib(ml_gen2);
            uint32_t x2 = distrib(ml_gen2);
            data.v_x_test.push_back( {x1, x2} );
            for( int j{0}; j<nBits; ++j )
            {
                if( ( x1 >> j ) & 1u == 1u )
                    kitty::set_bit( data.x_test[j], i );
                else
                    kitty::clear_bit( data.x_test[j], i );
                if( ( x2 >> j ) & 1u == 1u )
                    kitty::set_bit( data.x_test[j+nBits], i );
                else
                    kitty::clear_bit( data.x_test[j+nBits], i );
            }
            if( x2 > x1 )
                kitty::set_bit( data.y_test, i );
            else
                kitty::clear_bit( data.y_test, i );
            data.v_y_test.push_back( kitty::get_bit( data.y_test, i ) );
        }
        return data;
    }

    binary_classification_dataset_t linearly_separable_dataset_termometer( int nBits, int nTrain, int nValid, int nTest )
    {
        assert(pow(2,nBits) <= 32);
        binary_classification_dataset_t data;
        std::uniform_int_distribution<> distrib(0, ( 1u << nBits ) - 1 );
        uint32_t ORIGIN = 0xFFFFFFFF;
        PTT train0(nTrain);
        PTT valid0(nValid);
        PTT test0(nTest); 

        for( int i{0}; i<2*pow(2,nBits); ++i )
        {
            data.x_train.push_back( train0 );
            data.x_valid.push_back( valid0 );
            data.x_test.push_back( test0 );
        }

        data.y_train = train0;
        data.y_valid = valid0;
        data.y_test = test0;
        data.m_train = ~train0;
        data.m_valid = ~valid0;
        data.m_test  = ~test0;

        int i=0;
        std::set<std::pair<uint32_t, uint32_t>> set;
        int rep=0;
        while( i<nTrain && rep < 1000 )
        {
            rep++;
            uint32_t x1 = distrib(ml_gen2);
            uint32_t x2 = distrib(ml_gen2);
            uint32_t X1 = ORIGIN >> (32-x1);
            uint32_t X2 = ORIGIN >> (32-x2);

            std::pair<uint32_t, uint32_t> pair = std::make_pair( x1, x2 );
            if( set.find( pair ) != set.end() )
                continue;
            set.insert(pair);
            i++;
            data.v_x_train.push_back( {x1, x2} );
            for( int j{0}; j<pow(2,nBits); ++j )
            {
                if( ( X1 >> j ) & 1u == 1u )
                    kitty::set_bit( data.x_train[j], i );
                else
                    kitty::clear_bit( data.x_train[j], i );

                if( ( X2 >> j ) & 1u == 1u )
                    kitty::set_bit( data.x_train[j+nBits], i );
                else
                    kitty::clear_bit( data.x_train[j+nBits], i );
            }
            if( x2 > x1 )
                kitty::set_bit( data.y_train, i );
            else
                kitty::clear_bit( data.y_train, i );
            data.v_y_train.push_back( kitty::get_bit( data.y_train, i ) );
        }

        for( int i{0}; i<nValid; ++i )
        {
            uint32_t x1 = distrib(ml_gen2);
            uint32_t x2 = distrib(ml_gen2);
            uint32_t X1 = ORIGIN >> (32-x1);
            uint32_t X2 = ORIGIN >> (32-x2);

            data.v_x_valid.push_back( {x1, x2} );
            for( int j{0}; j<pow(2,nBits); ++j )
            {
                if( ( X1 >> j ) & 1u == 1u )
                    kitty::set_bit( data.x_valid[j], i );
                else
                    kitty::clear_bit( data.x_valid[j], i );
                if( ( X2 >> j ) & 1u == 1u )
                    kitty::set_bit( data.x_valid[j+nBits], i );
            }
            if( x2 > x1 )
                kitty::set_bit( data.y_valid, i );
            else
                kitty::clear_bit( data.y_valid, i );
            data.v_y_valid.push_back( kitty::get_bit( data.y_valid, i ) );
        }

        for( int i{0}; i<nTest; ++i )
        {
            uint32_t x1 = distrib(ml_gen2);
            uint32_t x2 = distrib(ml_gen2);
            uint32_t X1 = ORIGIN >> (32-x1);
            uint32_t X2 = ORIGIN >> (32-x2);
            
            data.v_x_test.push_back( {x1, x2} );
            for( int j{0}; j<pow(2,nBits); ++j )
            {
                if( ( X1 >> j ) & 1u == 1u )
                    kitty::set_bit( data.x_test[j], i );
                else
                    kitty::clear_bit( data.x_test[j], i );
                if( ( X2 >> j ) & 1u == 1u )
                    kitty::set_bit( data.x_test[j+nBits], i );
                else
                    kitty::clear_bit( data.x_test[j+nBits], i );
            }
            if( x2 > x1 )
                kitty::set_bit( data.y_test, i );
            else
                kitty::clear_bit( data.y_test, i );
            data.v_y_test.push_back( kitty::get_bit( data.y_test, i ) );
        }
        return data;
    }

} // namespace mcts

} // namespace mockturtle
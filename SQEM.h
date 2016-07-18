/*
 *     SQEM - Spherical Quadric Error Metric
 *
 *     Authors:  Jean-Marc THIERY, Emilie GUY and Tamy BOUBEKEUR
 *
 *     Copyright © 2012-2013 Telecom ParisTech - CNRS LTCI - Institut Mines Telecom
 *              All rights reserved
 *
 * This file is a part of the standalone implementation of the 
 * 
 * Spherical Quadric Error Metric.
 * 
 * For more information or if you use this file and need to reference it, please
 * 
 * refer to the following publication:
 *
 *    Sphere-Meshes: Shape Approximation using Spherical Quadric Error Metrics
 *    
 *    Jean-Marc Thiery, Emilie Guy and Tamy Boubekeur
 *    
 *    ACM Transaction on Graphics (Proc. SIGGRAPH Asia 2013), 32(6), Art. 178
 *
 *    http://www.telecom-paristech.fr/~boubek/papers/SphereMeshes/
 * 
 * SQEM is free software: you can redistribute it and/or modify
 *
 * it under the terms of the GNU Lesser General Public License as published by
 *
 * the Free Software Foundation, either version 3 of the License, or
 *
 * (at your option) any later version.
 *
 * SQEM is distributed in the hope that it will be useful,
 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef SQEM_H 
#define SQEM_H 

 #include <cmath>

static
unsigned int symmetric44IndicesToCompressedVectorIndex[4][4] =
{
    {0,1,2,3},
    {1,4,5,6},
    {2,5,7,8},
    {3,6,8,9}
};


// TODO:
// check for "dirty cases" where the cost is NaN
// (this should never happen in the wonderful world of mathematics, but it happens in the painful world of computer science)
// If you don't do it, and feed a priority_queue with NaN as a cost, it messes up completely the queue

// we minimize the quadratic energy 1/2 w^t * SQEM_A * w - SQEM_b * w:

/// The SQEM class is a 4D quadric which models the squared distance of a sphere to a set of oriented planes.
/// SQEMs can be summed up and a given SQEM can be used to optimize a sphere which approximate the set of planes.
class SQEM {
    // compressed version:
    double SQEM_A_compr[10];
    double SQEM_b_compr[4];
    double SQEM_c_compr;

    // the spherical quadric error is computed as x^t*A*x/2   -  b^t*x  +  c
    // A is stored as:
    // SQEM_A_compr[0]  SQEM_A_compr[1]  SQEM_A_compr[2]  SQEM_A_compr[3]
    // SQEM_A_compr[1]  SQEM_A_compr[4]  SQEM_A_compr[5]  SQEM_A_compr[6]
    // SQEM_A_compr[2]  SQEM_A_compr[5]  SQEM_A_compr[7]  SQEM_A_compr[8]
    // SQEM_A_compr[3]  SQEM_A_compr[6]  SQEM_A_compr[8]  SQEM_A_compr[9]

    SQEM( double a1 , double a2 , double a3 , double a4 , double a5 , double a6 , double a7 , double a8 , double a9 , double a10 ,
                 double b1 , double b2 , double b3 , double b4 ,
                 double c )
    {
        SQEM_A_compr[0] = a1;
        SQEM_A_compr[1] = a2;
        SQEM_A_compr[2] = a3;
        SQEM_A_compr[3] = a4;
        SQEM_A_compr[4] = a5;
        SQEM_A_compr[5] = a6;
        SQEM_A_compr[6] = a7;
        SQEM_A_compr[7] = a8;
        SQEM_A_compr[8] = a9;
        SQEM_A_compr[9] = a10;
        SQEM_b_compr[0] = b1;
        SQEM_b_compr[1] = b2;
        SQEM_b_compr[3] = b3;
        SQEM_b_compr[3] = b4;
        SQEM_c_compr = c;
    }

    inline
    double SQEM_A_determinant_with_column0_replaced( double x , double y , double z , double w ) const
    {
        return x * (
                    cSQEM_A(1,1) * ( cSQEM_A(2,2) * cSQEM_A(3,3) - cSQEM_A(3,2) * cSQEM_A(2,3) ) -
                    cSQEM_A(1,2) * ( cSQEM_A(2,1) * cSQEM_A(3,3) - cSQEM_A(3,1) * cSQEM_A(2,3) ) +
                    cSQEM_A(1,3) * ( cSQEM_A(2,1) * cSQEM_A(3,2) - cSQEM_A(3,1) * cSQEM_A(2,2) )
                    )  -
                cSQEM_A(0,1) * (
                    y * ( cSQEM_A(2,2) * cSQEM_A(3,3) - cSQEM_A(3,2) * cSQEM_A(2,3) ) -
                    cSQEM_A(1,2) * ( z * cSQEM_A(3,3) - w * cSQEM_A(2,3) ) +
                    cSQEM_A(1,3) * ( z * cSQEM_A(3,2) - w * cSQEM_A(2,2) )
                    )  +
                cSQEM_A(0,2) * (
                    y * ( cSQEM_A(2,1) * cSQEM_A(3,3) - cSQEM_A(3,1) * cSQEM_A(2,3) ) -
                    cSQEM_A(1,1) * ( z * cSQEM_A(3,3) - w * cSQEM_A(2,3) ) +
                    cSQEM_A(1,3) * ( z * cSQEM_A(3,1) - w * cSQEM_A(2,1) )
                    )  -
                cSQEM_A(0,3) * (
                    y * ( cSQEM_A(2,1) * cSQEM_A(3,2) - cSQEM_A(3,1) * cSQEM_A(2,2) ) -
                    cSQEM_A(1,1) * ( z * cSQEM_A(3,2) - w * cSQEM_A(2,2) ) +
                    cSQEM_A(1,2) * ( z * cSQEM_A(3,1) - w * cSQEM_A(2,1) )
                    );
    }

    inline
    double SQEM_A_determinant_with_column1_replaced( double x , double y , double z , double w ) const
    {
        return cSQEM_A(0,0) * (
                    y * ( cSQEM_A(2,2) * cSQEM_A(3,3) - cSQEM_A(3,2) * cSQEM_A(2,3) ) -
                    cSQEM_A(1,2) * ( z * cSQEM_A(3,3) - w * cSQEM_A(2,3) ) +
                    cSQEM_A(1,3) * ( z * cSQEM_A(3,2) - w * cSQEM_A(2,2) )
                    )  -
                x * (
                    cSQEM_A(1,0) * ( cSQEM_A(2,2) * cSQEM_A(3,3) - cSQEM_A(3,2) * cSQEM_A(2,3) ) -
                    cSQEM_A(1,2) * ( cSQEM_A(2,0) * cSQEM_A(3,3) - cSQEM_A(3,0) * cSQEM_A(2,3) ) +
                    cSQEM_A(1,3) * ( cSQEM_A(2,0) * cSQEM_A(3,2) - cSQEM_A(3,0) * cSQEM_A(2,2) )
                    )  +
                cSQEM_A(0,2) * (
                    cSQEM_A(1,0) * ( z * cSQEM_A(3,3) - w * cSQEM_A(2,3) ) -
                    y * ( cSQEM_A(2,0) * cSQEM_A(3,3) - cSQEM_A(3,0) * cSQEM_A(2,3) ) +
                    cSQEM_A(1,3) * ( cSQEM_A(2,0) * w - cSQEM_A(3,0) * z )
                    )  -
                cSQEM_A(0,3) * (
                    cSQEM_A(1,0) * ( z * cSQEM_A(3,2) - w * cSQEM_A(2,2) ) -
                    y * ( cSQEM_A(2,0) * cSQEM_A(3,2) - cSQEM_A(3,0) * cSQEM_A(2,2) ) +
                    cSQEM_A(1,2) * ( cSQEM_A(2,0) * w - cSQEM_A(3,0) * z )
                    );
    }

    inline
    double SQEM_A_determinant_with_column2_replaced( double x , double y , double z , double w ) const
    {
        return cSQEM_A(0,0) * (
                    cSQEM_A(1,1) * ( z * cSQEM_A(3,3) - w * cSQEM_A(2,3) ) -
                    y * ( cSQEM_A(2,1) * cSQEM_A(3,3) - cSQEM_A(3,1) * cSQEM_A(2,3) ) +
                    cSQEM_A(1,3) * ( cSQEM_A(2,1) * w - cSQEM_A(3,1) * z )
                    )  -
                cSQEM_A(0,1) * (
                    cSQEM_A(1,0) * ( z * cSQEM_A(3,3) - w * cSQEM_A(2,3) ) -
                    y * ( cSQEM_A(2,0) * cSQEM_A(3,3) - cSQEM_A(3,0) * cSQEM_A(2,3) ) +
                    cSQEM_A(1,3) * ( cSQEM_A(2,0) * w - cSQEM_A(3,0) * z )
                    )  +
                x * (
                    cSQEM_A(1,0) * ( cSQEM_A(2,1) * cSQEM_A(3,3) - cSQEM_A(3,1) * cSQEM_A(2,3) ) -
                    cSQEM_A(1,1) * ( cSQEM_A(2,0) * cSQEM_A(3,3) - cSQEM_A(3,0) * cSQEM_A(2,3) ) +
                    cSQEM_A(1,3) * ( cSQEM_A(2,0) * cSQEM_A(3,1) - cSQEM_A(3,0) * cSQEM_A(2,1) )
                    )  -
                cSQEM_A(0,3) * (
                    cSQEM_A(1,0) * ( cSQEM_A(2,1) * w - cSQEM_A(3,1) * z ) -
                    cSQEM_A(1,1) * ( cSQEM_A(2,0) * w - cSQEM_A(3,0) * z ) +
                    y * ( cSQEM_A(2,0) * cSQEM_A(3,1) - cSQEM_A(3,0) * cSQEM_A(2,1) )
                    );
    }

    inline
    double SQEM_A_determinant_with_column3_replaced( double x , double y , double z , double w ) const
    {
        return cSQEM_A(0,0) * (
                    cSQEM_A(1,1) * ( cSQEM_A(2,2) * w - cSQEM_A(3,2) * z ) -
                    cSQEM_A(1,2) * ( cSQEM_A(2,1) * w - cSQEM_A(3,1) * z ) +
                    y * ( cSQEM_A(2,1) * cSQEM_A(3,2) - cSQEM_A(3,1) * cSQEM_A(2,2) )
                    )  -
                cSQEM_A(0,1) * (
                    cSQEM_A(1,0) * ( cSQEM_A(2,2) * w - cSQEM_A(3,2) * z ) -
                    cSQEM_A(1,2) * ( cSQEM_A(2,0) * w - cSQEM_A(3,0) * z ) +
                    y * ( cSQEM_A(2,0) * cSQEM_A(3,2) - cSQEM_A(3,0) * cSQEM_A(2,2) )
                    )  +
                cSQEM_A(0,2) * (
                    cSQEM_A(1,0) * ( cSQEM_A(2,1) * w - cSQEM_A(3,1) * z ) -
                    cSQEM_A(1,1) * ( cSQEM_A(2,0) * w - cSQEM_A(3,0) * z ) +
                    y * ( cSQEM_A(2,0) * cSQEM_A(3,1) - cSQEM_A(3,0) * cSQEM_A(2,1) )
                    )  -
                x * (
                    cSQEM_A(1,0) * ( cSQEM_A(2,1) * cSQEM_A(3,2) - cSQEM_A(3,1) * cSQEM_A(2,2) ) -
                    cSQEM_A(1,1) * ( cSQEM_A(2,0) * cSQEM_A(3,2) - cSQEM_A(3,0) * cSQEM_A(2,2) ) +
                    cSQEM_A(1,2) * ( cSQEM_A(2,0) * cSQEM_A(3,1) - cSQEM_A(3,0) * cSQEM_A(2,1) )
                    );
    }

    inline
    double QEM_A_determinant_with_column0_replaced( double x , double y , double z ) const
    {
        return x * ( cSQEM_A(1,1) * cSQEM_A(2,2) - cSQEM_A(2,1) * cSQEM_A(1,2) ) -
                cSQEM_A(0,1) * ( y * cSQEM_A(2,2) - z * cSQEM_A(1,2) ) +
                cSQEM_A(0,2) * ( y * cSQEM_A(2,1) - z * cSQEM_A(1,1) );
    }


    inline
    double QEM_A_determinant_with_column1_replaced( double x , double y , double z ) const
    {
        return cSQEM_A(0,0) * ( y * cSQEM_A(2,2) - z * cSQEM_A(1,2) ) -
                x * ( cSQEM_A(1,0) * cSQEM_A(2,2) - cSQEM_A(2,0) * cSQEM_A(1,2) ) +
                cSQEM_A(0,2) * ( cSQEM_A(1,0) * z - cSQEM_A(2,0) * y );
    }

    inline
    double QEM_A_determinant_with_column2_replaced( double x , double y , double z ) const
    {
        return cSQEM_A(0,0) * ( cSQEM_A(1,1) * z - cSQEM_A(2,1) * y ) -
                cSQEM_A(0,1) * ( cSQEM_A(1,0) * z - cSQEM_A(2,0) * y ) +
                x * ( cSQEM_A(1,0) * cSQEM_A(2,1) - cSQEM_A(2,0) * cSQEM_A(1,1) );
    }


public:
    inline 
    SQEM() {}

    inline 
    void setZero()  {
        for( unsigned int i = 0 ; i < 10 ; ++i )
            SQEM_A_compr[i] = 0.0;
        for( unsigned int i = 0 ; i < 4 ; ++i )
            SQEM_b_compr[i] = 0.0;
        SQEM_c_compr = 0.0;
    }

    /// Necessary to compile but hopefully never used.
    /// Required since we use a std::priority_queue< std::pair< blabla , SQEM > >
    inline 
    bool operator < (const SQEM &  ) const { return true; }
    
    template< class point_t >
    inline 
    SQEM( const point_t & p , const point_t & n ) { setFromPlan(p,n); }

    template< class point_t >
    inline 
    void setFromPlan( const point_t & p , const point_t & n ) {
        double dot_product = p[0]*n[0] + p[1]*n[1] + p[2]*n[2];

        SQEM_A_compr[0] = 2.0*n[0]*n[0];
        SQEM_A_compr[1] = 2.0*n[0]*n[1];
        SQEM_A_compr[2] = 2.0*n[0]*n[2];
        SQEM_A_compr[3] = 2.0*n[0];
        SQEM_A_compr[4] = 2.0*n[1]*n[1];
        SQEM_A_compr[5] = 2.0*n[1]*n[2];
        SQEM_A_compr[6] = 2.0*n[1];
        SQEM_A_compr[7] = 2.0*n[2]*n[2];
        SQEM_A_compr[8] = 2.0*n[2];
        SQEM_A_compr[9] = 2.0;

        SQEM_b_compr[0] = 2.0 * dot_product * n[0];
        SQEM_b_compr[1] = 2.0 * dot_product * n[1];
        SQEM_b_compr[2] = 2.0 * dot_product * n[2];
        SQEM_b_compr[3] = 2.0 * dot_product;

        SQEM_c_compr = dot_product*dot_product;
    }

    inline 
    SQEM operator + (const SQEM & q2) const {
        SQEM res;

        for( unsigned int i = 0 ; i < 10 ; ++i )
            res.SQEM_A_compr[i] = this->SQEM_A_compr[i] + q2.SQEM_A_compr[i];
        for( unsigned int i = 0 ; i < 4 ; ++i )
            res.SQEM_b_compr[i] = this->SQEM_b_compr[i] + q2.SQEM_b_compr[i];
        res.SQEM_c_compr = this->SQEM_c_compr + q2.SQEM_c_compr;

        return res;
    }
    
    inline 
    SQEM operator * (double w) {
        SQEM res;

        for( unsigned int i = 0 ; i < 10 ; ++i )
            res.SQEM_A_compr[i] = this->SQEM_A_compr[i] * w;
        for( unsigned int i = 0 ; i < 4 ; ++i )
            res.SQEM_b_compr[i] = this->SQEM_b_compr[i] * w;
        res.SQEM_c_compr = this->SQEM_c_compr * w;

        return res;
    }
    
    inline 
    void operator *= (double w) {
        for( unsigned int i = 0 ; i < 10 ; ++i )
            this->SQEM_A_compr[i] = this->SQEM_A_compr[i] * w;
        for( unsigned int i = 0 ; i < 4 ; ++i )
            this->SQEM_b_compr[i] = this->SQEM_b_compr[i] * w;
        this->SQEM_c_compr = this->SQEM_c_compr * w;
    }
    
    inline 
    void operator += (const SQEM & q2) {
        for( unsigned int i = 0 ; i < 10 ; ++i )
            this->SQEM_A_compr[i] = this->SQEM_A_compr[i] + q2.SQEM_A_compr[i];
        for( unsigned int i = 0 ; i < 4 ; ++i )
            this->SQEM_b_compr[i] = this->SQEM_b_compr[i] + q2.SQEM_b_compr[i];
        this->SQEM_c_compr = this->SQEM_c_compr + q2.SQEM_c_compr;
    }

    inline
    double cSQEM_A( int i , int j ) const {
        return SQEM_A_compr[  symmetric44IndicesToCompressedVectorIndex[i][j]  ];
    }

    inline
    double cQEM_A( int i , int j ) const {
        return cSQEM_A(i,j);
    }

    inline
    double cSQEM_b( int i )const {
        return SQEM_b_compr[i];
    }

    inline
    double cQEM_b( int i ) const{
        return cSQEM_b(i);
    }

    inline
    double cSQEM_c( ) const{
        return SQEM_c_compr;
    }

    inline
    double cQEM_c( ) const{
        return cSQEM_c( );
    }

    template< class point_3D_t , class type_t >
    inline
    double evaluate( const point_3D_t & sphereCenter , type_t sphereRadius ) const {
        return 0.5 * (     cSQEM_A(0,0)*sphereCenter[0]*sphereCenter[0]+
                           2.0*cSQEM_A(0,1)*sphereCenter[0]*sphereCenter[1]+
                           2.0*cSQEM_A(0,2)*sphereCenter[0]*sphereCenter[2]+
                           2.0*cSQEM_A(0,3)*sphereCenter[0]*sphereRadius   +
                           cSQEM_A(1,1)*sphereCenter[1]*sphereCenter[1]+
                           2.0*cSQEM_A(1,2)*sphereCenter[1]*sphereCenter[2]+
                           2.0*cSQEM_A(1,3)*sphereCenter[1]*sphereRadius   +
                           cSQEM_A(2,2)*sphereCenter[2]*sphereCenter[2]+
                           2.0*cSQEM_A(2,3)*sphereCenter[2]*sphereRadius   +
                           cSQEM_A(3,3)*sphereRadius   *sphereRadius)
                - (cSQEM_b(0)*sphereCenter[0]+cSQEM_b(1)*sphereCenter[1]+cSQEM_b(2)*sphereCenter[2]+cSQEM_b(3)*sphereRadius)
                + cSQEM_c();
    }
    
    template< class point_4D_t >
    inline
    double evaluate( const point_4D_t & sphere ) const {
        return 0.5 * (     cSQEM_A(0,0)*sphere[0]*sphere[0]+
                           2.0*cSQEM_A(0,1)*sphere[0]*sphere[1]+
                           2.0*cSQEM_A(0,2)*sphere[0]*sphere[2]+
                           2.0*cSQEM_A(0,3)*sphere[0]*sphere[3]   +
                           cSQEM_A(1,1)*sphere[1]*sphere[1]+
                           2.0*cSQEM_A(1,2)*sphere[1]*sphere[2]+
                           2.0*cSQEM_A(1,3)*sphere[1]*sphere[3]   +
                           cSQEM_A(2,2)*sphere[2]*sphere[2]+
                           2.0*cSQEM_A(2,3)*sphere[2]*sphere[3]   +
                           cSQEM_A(3,3)*sphere[3]*sphere[3])
                - (cSQEM_b(0)*sphere[0]+cSQEM_b(1)*sphere[1]+cSQEM_b(2)*sphere[2]+cSQEM_b(3)*sphere[3])
                + cSQEM_c();
    }

    inline
    double SQEM_A_determinant() const
    {
        return cSQEM_A(0,0) * (
                    cSQEM_A(1,1) * ( cSQEM_A(2,2) * cSQEM_A(3,3) - cSQEM_A(3,2) * cSQEM_A(2,3) ) -
                    cSQEM_A(1,2) * ( cSQEM_A(2,1) * cSQEM_A(3,3) - cSQEM_A(3,1) * cSQEM_A(2,3) ) +
                    cSQEM_A(1,3) * ( cSQEM_A(2,1) * cSQEM_A(3,2) - cSQEM_A(3,1) * cSQEM_A(2,2) )
                    )  -
                cSQEM_A(0,1) * (
                    cSQEM_A(1,0) * ( cSQEM_A(2,2) * cSQEM_A(3,3) - cSQEM_A(3,2) * cSQEM_A(2,3) ) -
                    cSQEM_A(1,2) * ( cSQEM_A(2,0) * cSQEM_A(3,3) - cSQEM_A(3,0) * cSQEM_A(2,3) ) +
                    cSQEM_A(1,3) * ( cSQEM_A(2,0) * cSQEM_A(3,2) - cSQEM_A(3,0) * cSQEM_A(2,2) )
                    )  +
                cSQEM_A(0,2) * (
                    cSQEM_A(1,0) * ( cSQEM_A(2,1) * cSQEM_A(3,3) - cSQEM_A(3,1) * cSQEM_A(2,3) ) -
                    cSQEM_A(1,1) * ( cSQEM_A(2,0) * cSQEM_A(3,3) - cSQEM_A(3,0) * cSQEM_A(2,3) ) +
                    cSQEM_A(1,3) * ( cSQEM_A(2,0) * cSQEM_A(3,1) - cSQEM_A(3,0) * cSQEM_A(2,1) )
                    )  -
                cSQEM_A(0,3) * (
                    cSQEM_A(1,0) * ( cSQEM_A(2,1) * cSQEM_A(3,2) - cSQEM_A(3,1) * cSQEM_A(2,2) ) -
                    cSQEM_A(1,1) * ( cSQEM_A(2,0) * cSQEM_A(3,2) - cSQEM_A(3,0) * cSQEM_A(2,2) ) +
                    cSQEM_A(1,2) * ( cSQEM_A(2,0) * cSQEM_A(3,1) - cSQEM_A(3,0) * cSQEM_A(2,1) )
                    );
    }
    
    inline
    double QEM_A_determinant() const
    {
        return cSQEM_A(0,0) * ( cSQEM_A(1,1) * cSQEM_A(2,2) - cSQEM_A(2,1) * cSQEM_A(1,2) ) -
                cSQEM_A(0,1) * ( cSQEM_A(1,0) * cSQEM_A(2,2) - cSQEM_A(2,0) * cSQEM_A(1,2) ) +
                cSQEM_A(0,2) * ( cSQEM_A(1,0) * cSQEM_A(2,1) - cSQEM_A(2,0) * cSQEM_A(1,1) );
    }

    /// This function minimizes the quadric, with no maximum radius constraint,
    /// and ensures that the center is on the segment [pa,pb] in case it is degenerate.
    template< class point_t , class type_t >
    double minimize( point_t & sphereCenter , type_t & sphereRadius , const point_t & pa , const point_t & pb ) const {
        double det_epsilon = 0.0001;
        // 1) check if SQEM_A is invertible:
        double det44 = SQEM_A_determinant();
        if( std::abs(det44) > det_epsilon ) {
            sphereRadius = SQEM_A_determinant_with_column3_replaced( cSQEM_b(0) , cSQEM_b(1) , cSQEM_b(2) , cSQEM_b(3) ) / det44;

            // 2) check if it is a convex sphere and is in the half space that is allowed:
            if( sphereRadius >= 0.0 ) {
                sphereCenter[0] = SQEM_A_determinant_with_column0_replaced( cSQEM_b(0) , cSQEM_b(1) , cSQEM_b(2) , cSQEM_b(3) ) / det44;
                sphereCenter[1] = SQEM_A_determinant_with_column1_replaced( cSQEM_b(0) , cSQEM_b(1) , cSQEM_b(2) , cSQEM_b(3) ) / det44;
                sphereCenter[2] = SQEM_A_determinant_with_column2_replaced( cSQEM_b(0) , cSQEM_b(1) , cSQEM_b(2) , cSQEM_b(3) ) / det44;
                return evaluate( sphereCenter , sphereRadius );
            }
        }

        // if sphereRadius < 0 for the global minimizer, then the minimizer on the restriction r >= 0 is found on the hyperplane r = 0:
        sphereRadius = 0.0;
        double det33 = QEM_A_determinant();
        if( std::abs(det33) > det_epsilon ) {
            sphereCenter[0] = QEM_A_determinant_with_column0_replaced(cSQEM_b(0) , cSQEM_b(1) , cSQEM_b(2)) / det33;
            sphereCenter[1] = QEM_A_determinant_with_column1_replaced(cSQEM_b(0) , cSQEM_b(1) , cSQEM_b(2)) / det33;
            sphereCenter[2] = QEM_A_determinant_with_column2_replaced(cSQEM_b(0) , cSQEM_b(1) , cSQEM_b(2)) / det33;
            return evaluate(sphereCenter , sphereRadius);
        }

        // - otherwise evaluate Q(x) on the segment:
        {
            const point_t & ab = pb-pa;
            point_t Aab(
                        cSQEM_A(0,0) * ab[0] + cSQEM_A(0,1) * ab[1] + cSQEM_A(0,2) * ab[2] ,
                        cSQEM_A(1,0) * ab[0] + cSQEM_A(1,1) * ab[1] + cSQEM_A(1,2) * ab[2] ,
                        cSQEM_A(2,0) * ab[0] + cSQEM_A(2,1) * ab[1] + cSQEM_A(2,2) * ab[2]);
            double abAab = Aab[0] * ab[0] + Aab[1] * ab[1] + Aab[2] * ab[2];
            double lambda_minimizer = ((ab[0] * pa[0] + ab[1] * pa[1] + ab[2] * pa[2])   -   (Aab[0] * pa[0] + Aab[1] * pa[1] + Aab[2] * pa[2])) / abAab;
            if( lambda_minimizer >= 0.0  &&  lambda_minimizer <= 1.0 ) {
                sphereCenter = pa + lambda_minimizer * (pb - pa);
                return evaluate(sphereCenter , sphereRadius);
            }
        }

        // - if it fails then return the midpoint.
        sphereCenter = (pb + pa) / 2.f;
        return evaluate(sphereCenter , sphereRadius);
    }

    /// This function takes into account basic linear inequalities or the form r >= 0  &&  r <= MAX_RADIUS,
    /// and ensures that the center is on the segment [pa,pb] in case it is degenerate
    /// (if we really don't know what to do, we put it at the center of the segment):
    template< class point_t , class type_t >
    double minimize( point_t & sphereCenter , type_t & sphereRadius , const point_t & pa , const point_t & pb , double MAX_RADIUS ) const {
        // we minimize the quadratic energy 1/2 w^t * SQEM_A * w - SQEM_b * w:
        double det_epsilon = 0.0001;
        // 1) check if SQEM_A is invertible:
        double det44 = SQEM_A_determinant();
        if( std::abs(det44) > det_epsilon ) {
            sphereRadius = SQEM_A_determinant_with_column3_replaced( cSQEM_b(0) , cSQEM_b(1) , cSQEM_b(2) , cSQEM_b(3) ) / det44;

            // 2) check if it is a convex sphere and is in the half space that is allowed:
            if( sphereRadius >= 0.0   &&   sphereRadius  <= MAX_RADIUS ) {
                sphereCenter[0] = SQEM_A_determinant_with_column0_replaced( cSQEM_b(0) , cSQEM_b(1) , cSQEM_b(2) , cSQEM_b(3) ) / det44;
                sphereCenter[1] = SQEM_A_determinant_with_column1_replaced( cSQEM_b(0) , cSQEM_b(1) , cSQEM_b(2) , cSQEM_b(3) ) / det44;
                sphereCenter[2] = SQEM_A_determinant_with_column2_replaced( cSQEM_b(0) , cSQEM_b(1) , cSQEM_b(2) , cSQEM_b(3) ) / det44;
                return evaluate( sphereCenter , sphereRadius );
            }
            else
            {
                // 3) if not, then the minimizer is either on the hyperplan r = 0 or r = MAX_RADIUS, compare them:
                double det33 = QEM_A_determinant();
                if( std::abs(det33) > det_epsilon ) {
                    double scoreOnNullRadiusHyperplan = 0.0;
                    double scoreOnMaxRadiusHyperplan = 0.0;

                    point_t minimizerCenterOnNullRadiusHyperplan, minimizerCenterOnMaxRadiusHyperplan;

                    point_t QEM_b , QEM_b_biased;
                    for( unsigned int i = 0 ; i < 3 ; ++i )
                        QEM_b[i] = cSQEM_b(i);
                    QEM_b_biased[0] = - SQEM_A_compr[3];
                    QEM_b_biased[1] = - SQEM_A_compr[6];
                    QEM_b_biased[2] = - SQEM_A_compr[8];

                    {
                        //                    minimizerCenterOnNullRadiusHyperplan = QEM_A.ldlt().solve(QEM_b);
                        minimizerCenterOnNullRadiusHyperplan[0] = QEM_A_determinant_with_column0_replaced( QEM_b[0] , QEM_b[1] , QEM_b[2] ) / det33;
                        minimizerCenterOnNullRadiusHyperplan[1] = QEM_A_determinant_with_column1_replaced( QEM_b[0] , QEM_b[1] , QEM_b[2] ) / det33;
                        minimizerCenterOnNullRadiusHyperplan[2] = QEM_A_determinant_with_column2_replaced( QEM_b[0] , QEM_b[1] , QEM_b[2] ) / det33;
                        scoreOnNullRadiusHyperplan = evaluate(minimizerCenterOnNullRadiusHyperplan , 0);
                    }
                    {
                        //                    minimizerCenterOnMaxRadiusHyperplan = QEM_A.ldlt().solve(QEM_b + MAX_RADIUS * QEM_b_biased);
                        minimizerCenterOnMaxRadiusHyperplan[0] = QEM_A_determinant_with_column0_replaced( QEM_b[0] + MAX_RADIUS * QEM_b_biased[0], QEM_b[1] + MAX_RADIUS * QEM_b_biased[1] , QEM_b[2] + MAX_RADIUS * QEM_b_biased[2] ) / det33;
                        minimizerCenterOnMaxRadiusHyperplan[1] = QEM_A_determinant_with_column1_replaced( QEM_b[0] + MAX_RADIUS * QEM_b_biased[0], QEM_b[1] + MAX_RADIUS * QEM_b_biased[1] , QEM_b[2] + MAX_RADIUS * QEM_b_biased[2] ) / det33;
                        minimizerCenterOnMaxRadiusHyperplan[2] = QEM_A_determinant_with_column2_replaced( QEM_b[0] + MAX_RADIUS * QEM_b_biased[0], QEM_b[1] + MAX_RADIUS * QEM_b_biased[1] , QEM_b[2] + MAX_RADIUS * QEM_b_biased[2] ) / det33;
                        scoreOnMaxRadiusHyperplan = evaluate(minimizerCenterOnMaxRadiusHyperplan , MAX_RADIUS);
                    }

                    if( scoreOnNullRadiusHyperplan  <=  scoreOnMaxRadiusHyperplan )
                    {
                        sphereRadius = 0.0;
                        sphereCenter[0] = minimizerCenterOnNullRadiusHyperplan[0];
                        sphereCenter[1] = minimizerCenterOnNullRadiusHyperplan[1];
                        sphereCenter[2] = minimizerCenterOnNullRadiusHyperplan[2];
                        return scoreOnNullRadiusHyperplan;
                    }
                    else
                    {
                        sphereRadius = MAX_RADIUS;
                        sphereCenter[0] = minimizerCenterOnMaxRadiusHyperplan[0];
                        sphereCenter[1] = minimizerCenterOnMaxRadiusHyperplan[1];
                        sphereCenter[2] = minimizerCenterOnMaxRadiusHyperplan[2];
                        return scoreOnMaxRadiusHyperplan;
                    }
                }
            }
        }


        /// If SQEM_A is NOT invertible, then we consider a sphere on the segment [va vb]:
        /// i) check if QEM_A is invertible:
        const point_t & ab = pb - pa;

        double A22_00 =
                (cSQEM_A(0,0) * ab[0] + cSQEM_A(0,1) * ab[1] + cSQEM_A(0,2) * ab[2])*ab[0] +
                (cSQEM_A(1,0) * ab[0] + cSQEM_A(1,1) * ab[1] + cSQEM_A(1,2) * ab[2])*ab[1] +
                (cSQEM_A(2,0) * ab[0] + cSQEM_A(2,1) * ab[1] + cSQEM_A(2,2) * ab[2])*ab[2]  ,
                A22_11 = cSQEM_A(3,3)   ,
                A22_01 = cSQEM_A(3,0) * ab[0] + cSQEM_A(3,1) * ab[1] + cSQEM_A(3,2) * ab[2];
        double A22_10 = A22_01;


        double b2_0 = cSQEM_b(0) * ab[0] + cSQEM_b(1) * ab[1] + cSQEM_b(2) * ab[2]
                - ( (cSQEM_A(0,0) * ab[0] + cSQEM_A(0,1) * ab[1] + cSQEM_A(0,2) * ab[2])*pa[0] +
                    (cSQEM_A(1,0) * ab[0] + cSQEM_A(1,1) * ab[1] + cSQEM_A(1,2) * ab[2])*pa[1] +
                    (cSQEM_A(2,0) * ab[0] + cSQEM_A(2,1) * ab[1] + cSQEM_A(2,2) * ab[2])*pa[2] );
        double b2_1 = cSQEM_b(3)
                - (cSQEM_A(3,0) * pa[0] + cSQEM_A(3,1) * pa[1] + cSQEM_A(3,2) * pa[2]);

        double c2 = cQEM_c() - (cSQEM_b(0) * pa[0] + cSQEM_b(1) * pa[1] + cSQEM_b(2) * pa[2]) + 0.5 * ( (cSQEM_A(0,0) * pa[0] + cSQEM_A(0,1) * pa[1] + cSQEM_A(0,2) * pa[2])*pa[0] +
                                                                                                        (cSQEM_A(1,0) * pa[0] + cSQEM_A(1,1) * pa[1] + cSQEM_A(1,2) * pa[2])*pa[1] +
                                                                                                        (cSQEM_A(2,0) * pa[0] + cSQEM_A(2,1) * pa[1] + cSQEM_A(2,2) * pa[2])*pa[2]) ;

        double det22 = A22_00 * A22_11  -  A22_01 * A22_10;
        if( std::abs(det22) > det_epsilon ) {
            double res_0 = (b2_0 * A22_11  -  b2_1 * A22_10) / det22;
            double res_1 = (A22_00 * b2_1  -  A22_01 * b2_0) / det22;

            if( res_0 >= 0.f   &&   res_0 <= 1.f    &&    res_1 >= 0.f   &&   res_1 <= MAX_RADIUS ) {
                sphereCenter = pa + res_0 * (pb - pa);
                sphereRadius = res_1;
                return evaluate(sphereCenter , sphereRadius);
            }

            // otherwise it means that the minimizer is either on lambda = 0 , lambda = 1 , r = 0 , or r = MAX_RADIUS:

            // we fix the position, and we find r_optimal between 0 and MAX_RADIUS:
            double minCost;
            // i: lambda = 0:
            {
                double lambda = 0;
                double rr = (b2_1 - A22_10 * lambda) / A22_11;
                if( rr <= 0.0 )
                {
                    rr = 0.0;
                }
                else if( rr >= MAX_RADIUS )
                {
                    rr = MAX_RADIUS;
                }

                double ccost = 0.5 * (  A22_00 * lambda * lambda   +   2.0 * lambda * rr * A22_10    +   A22_11 * rr * rr )
                        - ( b2_0 * lambda   +  b2_1 * rr )
                        + c2;

                sphereRadius = rr;
                sphereCenter = pa + lambda * (pb - pa);

                minCost = ccost;
            }
            // ii: lambda = 1:
            {
                double lambda = 1;
                double rr = (b2_1 - A22_10 * lambda) / A22_11;
                if( rr <= 0.0 )
                {
                    rr = 0.0;
                }
                else if( rr >= MAX_RADIUS )
                {
                    rr = MAX_RADIUS;
                }

                double ccost = 0.5 * (  A22_00 * lambda * lambda   +   2.0 * lambda * rr * A22_10    +   A22_11 * rr * rr )
                        - ( b2_0 * lambda   +  b2_1 * rr )
                        + c2;

                if( ccost < minCost ) {
                    sphereRadius = rr;
                    sphereCenter = pa + lambda * (pb - pa);

                    minCost = ccost;
                }
            }

            // we fix the radius, and we find the optimal position along ab:
            // iii: r = 0:
            {
                double rr = 0.0;
                double lambda = (b2_0 - A22_10 * rr) / A22_00;
                if( lambda <= 0.0 )
                {
                    lambda = 0.0;
                }
                else if( lambda >= 1.0 )
                {
                    lambda = 1.0;
                }

                double ccost = 0.5 * (  A22_00 * lambda * lambda   +   2.0 * lambda * rr * A22_10    +   A22_11 * rr * rr )
                        - ( b2_0 * lambda   +  b2_1 * rr )
                        + c2;

                if( ccost < minCost ) {
                    sphereRadius = rr;
                    sphereCenter = pa + lambda * (pb - pa);

                    minCost = ccost;
                }
            }
            // iv: r = MAX_RADIUS:
            {
                double rr = MAX_RADIUS;
                double lambda = (b2_0 - A22_10 * rr) / A22_00;
                if( lambda <= 0.0 )
                {
                    lambda = 0.0;
                }
                else if( lambda >= 1.0 )
                {
                    lambda = 1.0;
                }

                double ccost = 0.5 * (  A22_00 * lambda * lambda   +   2.0 * lambda * rr * A22_10    +   A22_11 * rr * rr )
                        - ( b2_0 * lambda   +  b2_1 * rr )
                        + c2;

                if( ccost < minCost ) {
                    sphereRadius = rr;
                    sphereCenter = pa + lambda * (pb - pa);

                    minCost = ccost;
                }
            }

            return minCost;
        }


        // If everything fails then return the midpoint (lambda = 0.5) and optimize the radius in [0 ; MAX_RADIUS]:
        // lambda = 0.5:
        {
            double lambda = 0.5;
            double rr = (b2_1 - A22_10 * lambda) / A22_11;
            if( rr <= 0.0 )
            {
                rr = 0.0;
            }
            else if( rr >= MAX_RADIUS )
            {
                rr = MAX_RADIUS;
            }

            double ccost = 0.5 * (  A22_00 * lambda * lambda   +   2.0 * lambda * rr * A22_10    +   A22_11 * rr * rr )
                    - ( b2_0 * lambda   +  b2_1 * rr )
                    + c2;

            sphereRadius = rr;
            sphereCenter = (pb + pa)/2.0;

            return ccost;
        }

        // get rid of warning:
        return 0.0;
    }

    /// This function takes into account basic linear equality or the form r = sphereRadius,
    /// and ensures that the center is on the segment [pa,pb] in case it is degenerate
    /// (if we really don't know what to do, we put it at the center of the segment):
    template< class point_t , class type_t >
    double minimizeWithFixedRadius( point_t & sphereCenter , type_t sphereRadius , const point_t & pa , const point_t & pb ) const {
        // we minimize the quadratic energy 1/2 w^t * SQEM_A * w - SQEM_b * w:

        double det_epsilon = 0.0001;

        // 1) check if SQEM_A is invertible:
        double det33 = QEM_A_determinant();
        if( std::abs(det33) > det_epsilon ) {
            double scoreOnFixedRadiusHyperplan = 0.0;

            point_t QEM_b , QEM_b_biased;
            for( unsigned int i = 0 ; i < 3 ; ++i )
                QEM_b[i] = cSQEM_b(i);
            QEM_b_biased[0] = - SQEM_A_compr[3];
            QEM_b_biased[1] = - SQEM_A_compr[6];
            QEM_b_biased[2] = - SQEM_A_compr[8];

            {
                sphereCenter[0] = QEM_A_determinant_with_column0_replaced( QEM_b[0] + sphereRadius * QEM_b_biased[0], QEM_b[1] + sphereRadius * QEM_b_biased[1] , QEM_b[2] + sphereRadius * QEM_b_biased[2] ) / det33;
                sphereCenter[1] = QEM_A_determinant_with_column1_replaced( QEM_b[0] + sphereRadius * QEM_b_biased[0], QEM_b[1] + sphereRadius * QEM_b_biased[1] , QEM_b[2] + sphereRadius * QEM_b_biased[2] ) / det33;
                sphereCenter[2] = QEM_A_determinant_with_column2_replaced( QEM_b[0] + sphereRadius * QEM_b_biased[0], QEM_b[1] + sphereRadius * QEM_b_biased[1] , QEM_b[2] + sphereRadius * QEM_b_biased[2] ) / det33;
                scoreOnFixedRadiusHyperplan = evaluate(sphereCenter , sphereRadius);
                return scoreOnFixedRadiusHyperplan;
            }
        }


        // If SQEM_A is NOT invertible, then we consider a sphere on the segment [va vb]:
        // i) check if QEM_A is invertible:
        const point_t & ab = pb - pa;

        double A22_00 =
                (cSQEM_A(0,0) * ab[0] + cSQEM_A(0,1) * ab[1] + cSQEM_A(0,2) * ab[2])*ab[0] +
                (cSQEM_A(1,0) * ab[0] + cSQEM_A(1,1) * ab[1] + cSQEM_A(1,2) * ab[2])*ab[1] +
                (cSQEM_A(2,0) * ab[0] + cSQEM_A(2,1) * ab[1] + cSQEM_A(2,2) * ab[2])*ab[2]  ,
                A22_11 = cSQEM_A(3,3)   ,
                A22_01 = cSQEM_A(3,0) * ab[0] + cSQEM_A(3,1) * ab[1] + cSQEM_A(3,2) * ab[2];
        double A22_10 = A22_01;


        double b2_0 = cSQEM_b(0) * ab[0] + cSQEM_b(1) * ab[1] + cSQEM_b(2) * ab[2]
                - ( (cSQEM_A(0,0) * ab[0] + cSQEM_A(0,1) * ab[1] + cSQEM_A(0,2) * ab[2])*pa[0] +
                    (cSQEM_A(1,0) * ab[0] + cSQEM_A(1,1) * ab[1] + cSQEM_A(1,2) * ab[2])*pa[1] +
                    (cSQEM_A(2,0) * ab[0] + cSQEM_A(2,1) * ab[1] + cSQEM_A(2,2) * ab[2])*pa[2] );
        double b2_1 = cSQEM_b(3)
                - (cSQEM_A(3,0) * pa[0] + cSQEM_A(3,1) * pa[1] + cSQEM_A(3,2) * pa[2]);

        double c2 = cQEM_c() - (cSQEM_b(0) * pa[0] + cSQEM_b(1) * pa[1] + cSQEM_b(2) * pa[2]) + 0.5 * ( (cSQEM_A(0,0) * pa[0] + cSQEM_A(0,1) * pa[1] + cSQEM_A(0,2) * pa[2])*pa[0] +
                                                                                                        (cSQEM_A(1,0) * pa[0] + cSQEM_A(1,1) * pa[1] + cSQEM_A(1,2) * pa[2])*pa[1] +
                                                                                                        (cSQEM_A(2,0) * pa[0] + cSQEM_A(2,1) * pa[1] + cSQEM_A(2,2) * pa[2])*pa[2]) ;

        double det22 = A22_00 * A22_11  -  A22_01 * A22_10;
        if( std::abs(det22) > det_epsilon ) {
            // otherwise it means that the minimizer is either on lambda = 0 , lambda = 1 , with r = sphereRadius:
            // we fix the radius, and we find lambda_optimal between 0 and 1:

            double minCost;
            // r = sphereRadius:
            {
                double rr = sphereRadius;
                double lambda = (b2_0 - A22_10 * rr) / A22_00;
                if( lambda <= 0.0 )
                {
                    lambda = 0.0;
                }
                else if( lambda >= 1.0 )
                {
                    lambda = 1.0;
                }
                // this simple test is valid, if lambdaMinimizer < 0, we know for sure that f(0) <= f(1) since it's a positive quadric

                double ccost = 0.5 * (  A22_00 * lambda * lambda   +   2.0 * lambda * rr * A22_10    +   A22_11 * rr * rr )
                        - ( b2_0 * lambda   +  b2_1 * rr )
                        + c2;

                {
                    sphereCenter = pa + lambda * (pb - pa);
                    minCost = ccost;
                }
            }

            return minCost;
        }

        // If it fails then return the midpoint (lambda = 0.5) and optimize the radius in [0 ; MAX_RADIUS]:
        // lambda = 0.5:
        {
            double lambda = 0.5;
            double rr = sphereRadius;
            double ccost = 0.5 * (  A22_00 * lambda * lambda   +   2.0 * lambda * rr * A22_10    +   A22_11 * rr * rr )
                    - ( b2_0 * lambda   +  b2_1 * rr )
                    + c2;

            sphereCenter = (pb + pa)/2.0;
            return ccost;
        }

        // get rid of warning:
        return 0.0;
    }
};



inline std::ostream & operator << (std::ostream & s , SQEM const & p)
{
    s << "SQEM_A:" << std::endl;
    s << p.cSQEM_A(0,0) << " \t" <<  p.cSQEM_A(0,1) << " \t" << p.cSQEM_A(0,2) << " \t" <<  p.cSQEM_A(0,3) << std::endl;
    s << p.cSQEM_A(1,0) << " \t" <<  p.cSQEM_A(1,1) << " \t" << p.cSQEM_A(1,2) << " \t" <<  p.cSQEM_A(1,3) << std::endl;
    s << p.cSQEM_A(2,0) << " \t" <<  p.cSQEM_A(2,1) << " \t" << p.cSQEM_A(2,2) << " \t" <<  p.cSQEM_A(2,3) << std::endl;
    s << p.cSQEM_A(3,0) << " \t" <<  p.cSQEM_A(3,1) << " \t" << p.cSQEM_A(3,2) << " \t" <<  p.cSQEM_A(3,3) << std::endl;
    s << "QEM_A:" << std::endl;
    s << p.cQEM_A(0,0) << " \t" <<  p.cQEM_A(0,1) << " \t" << p.cQEM_A(0,2) <<  std::endl;
    s << p.cQEM_A(1,0) << " \t" <<  p.cQEM_A(1,1) << " \t" << p.cQEM_A(1,2) <<  std::endl;
    s << p.cQEM_A(2,0) << " \t" <<  p.cQEM_A(2,1) << " \t" << p.cQEM_A(2,2) <<  std::endl;
    s << "SQEM_b:" << std::endl;
    s << p.cSQEM_b(0) << " \t" <<  p.cSQEM_b(1) << " \t" <<  p.cSQEM_b(2) << " \t" <<  p.cSQEM_b(3) << std::endl;
    s << "QEM_b:" << std::endl;
    s << p.cQEM_b(0) << " \t" <<  p.cQEM_b(1) << " \t" <<  p.cQEM_b(2) << std::endl;
    s << "QEM_c:" << std::endl;
    s << p.cSQEM_c() << std::endl << std::endl;
    return s;
}

#endif // SQEM_H 
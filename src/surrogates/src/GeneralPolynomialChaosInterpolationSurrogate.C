// This class
#include <queso/GeneralPolynomialChaosInterpolationSurrogate.h>
//Queso
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/MultiDimensionalIndexing.h>
#include <queso/InterpolationSurrogateData.h>

#include<stdio.h>
#include<gls/gsl_poly.h>

namespace QUESO
{


 template<class V, class M>
  GeneralPolynomialChaosInterpolationSurrogate<V,M>::GeneralPolynomialChaosInterpolationSurrogate(const InterpolationSurrogateData<V,M>& data)
    : InterpolationSurrogateBase<V,M>(data)
  {}

  template<class V, class M>
  double GeneralPolynomialChaosInterpolationSurrogate<V,M>::evaluate(const V & domainVector, unsigned int & dim) const
  {

    // lower bound coordinates at the current dimemsion
    std::vector<double> x_min(this->m_data.dim());

    // upper bound coordinates at the current dimension
    std::vector<double> x_max(this->m_data.dim());

    //sanity check
    queso_assert_equal_to( order.dim(), m_data.dim()):



    for( unsigned int d = 0; d < this->m_data.dim(); d++ )
      {

        x_min[d] = this->m_data.x_min(d);
        x_max[d] = this->m_data.x_max(d);
      
      }

    max_degree = order(dim);//this is to which order we want to build the hermite polynomial at the current dimension

 		//this is the matrix for the hermite polynimals. row number is the order, columnun is the coefficients for x^j
      	//Hi = hermite_coef(i)(j)*x^j for j = 0 to max_degree
  	    //for example
  	    //H0 =  1				hermite coef (2)(2) = [1  ,0  ,0  ]
  	    //H1 =  0 + 2x								  [0  ,2  ,0  ]
  	    //H2 = -2 + 0x + 4x^2,  					  [-2 ,0  ,4  ]
      	//initialize here
 	int hermite_coef[max_degree+1][max_degree+1] = zeros[max_degree+1][max_degree+1];

 	for (i = 0; i < max_degree+1; i++)
 	{
 		if (i == 0)
 		{
 			hermite_coef[0][0] = 1;//H0 = 1
 		}
 		else if (i == 1)
 		{
 			hermite_coef[0][0] = 0;
 			hermite_coef[1][1] = 2; //H1 = 0 + 2x
 		}
 		else
 		{
 			for (j = max_degree; j >= 0; j--)//recurrence law for hermites H_{n+1}(x) = 2x*H_n(x)-2n*H_{n-1}(x)
 			{
 				hermite_coef[i][j] = hermite_coef(i-1)(j-1)*2
 				hermite_coef[i][j] -= (i-1)hermite_coef(i-2)(j) 
 			}
 		}
 	}//finished constructing the hermite polynomial matrix

 	int n_quad_point = (max_degree+2)(max_degree-1)/2; //number of quadrature points in total is 2 + 3 +...+ n 

 	std::vector<double> qudrature_points_gaussian_hermite(n_quad_point);

 	//calculating quadrature points, from 2nd degree
 	for (unsigned int i = 2; i < max_degree + 1; i++)
 	{
 		int hermite[i+1] = zeros[i+1];//initializing the hermite poly nomial that we are currently working on 
 		for (unsigned int j = 0; j < max_degree + 1; j++)
 		{
 			if(hermite_coef[i][j] != 0)
 			{
 			hermite[j] = hermite_coef[i][j]; //polulating the hermite poly, and omiting the zeros at the end to ensre the leading term is non zeor
 			}
 		}
 		//now we have a array of the coefficients of the hermite poly of the current order that we are working on, staring from 0th term, then x term, x^2 term...
 		
 		//actual calculation starts, reference page:  https://www.gnu.org/software/gsl/manual/html_node/Roots-of-Polynomials-Examples.html#Roots-of-Polynomials-Examples
 		double q [2*i];
 		gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (i+1);
 		gsl_poly_complex_solve (a, i+1, w, q);
 		gsl_poly_complex_workspace_free (w);
 		//q[2(i)] is the roots for the hermite poly, alternating between real and imaginary, we can do a sanity check
 			for (unsigned int i = 0, i < sizeof(q); i++)
			{
				if (i%2 != 0 && q(i) != 0)
				{
					 queso_error	/* should have all real roots */
				}
			}
		//what we want are the quadrature points, they are on q[0], q[2], q[4]...
		//now we need to populate the values to our quadrature point array
		
			for (unsigned int b = 0; b < sizeof(q)/2; b++)
			{
				qudrature_points_gaussian_hermite[(2+i-1)*(i-2)/2+b] = q[b*2]
			}

		//now this is a bit messy, illustration is as follows:
		//
		//	i = 2			q[0]				q[2]
		//					 |		 	 		 |
		//					 V       	 		 V
		//				 q_p_g_h[0]		     q_p_g_h[1]
		//
		//
		//
		//	i = 3			q[0]				q[2]          q[4]
		//					 |		 	 		 |				|
		//					 V       	 		 V 				V
		//				 q_p_g_h[2]		     q_p_g_h[3]      q_p_g_h[4]
	 	//
		//
		//
	 	//	i = 4			q[0]				q[2]          q[4]				q[6]
		//					 |		 	 		 |				|				 |
		//					 V       	 		 V 				V 				 V
		//				 q_p_g_h[5]		     q_p_g_h[6]      q_p_g_h[7]		 q_p_g_h[8]
 	


 	}//now we have : hermite_coef[][], qudrature_points_gaussian_hermite[]




    // Now evaluate the interpolant
    return this->eval_interpolant( x_min, x_max, values, domainVector );
  }

 // template<class V, class M>
 //  void GeneralPolynomialChaosInterpolationSurrogate<V,M>::compute_interval_indices(const V & domainVector,
 //                                                                           std::vector<unsigned int>& indices) const
 //  {
 //    queso_assert_equal_to( domainVector.sizeGlobal(), this->m_data.dim() );
 //    queso_assert_equal_to( indices.size(), this->m_data.dim() );

 //    for( unsigned int d = 0; d < this->m_data.dim(); d++ )
 //      {
 //        double spacing = this->m_data.spacing(d);
 //        indices[d] = std::floor( (domainVector[d] - this->m_data.x_min(d))/spacing );

 //        // Index should be less than the number of point along this dimension
 //        queso_assert_less( indices[d], this->m_data.get_n_points()[d] );
 //      }
 //  }	

template<class V, class M>
   GeneralPolynomialChaosInterpolationSurrogate<V,M>::compute_quadrature_points(const std::vector<int>& order;) const
  {

  }


// template<class V, class M>
//   void GeneralPolynomialChaosInterpolationSurrogate<V,M>::compute_interval_values( const std::vector<unsigned int>& indices,
//                                                                            std::vector<double>& x_min,
//                                                                            std::vector<double>& x_max,
//                                                                            std::vector<double>& values ) const
//   {
//     queso_assert_equal_to( x_min.size(), this->m_data.dim() );
//     queso_assert_equal_to( x_max.size(), this->m_data.dim() );
//     // queso_assert_equal_to( values.size(), this->n_coeffs() );
//     // queso_assert_equal_to( indices.size(), this->m_data.dim() );

//     // // First, use the lower bound (global) indices to populate x_min, x_max
//     // for( unsigned int d = 0; d < this->m_data.dim(); d++ )
//     //   {
//     //     x_min[d] = this->m_data.get_x( d, indices[d] );
//     //     x_max[d] = x_min[d] + this->m_data.spacing( d );
//     //   }

//     // Now populate the values.
//     std::vector<unsigned int> local_indices(this->m_data.dim());
//     std::vector<unsigned int> global_indices(this->m_data.dim());
//     for( unsigned int n = 0 ; n < this->n_coeffs(); n++ )
//       {
//         // Figure out local indices (each coordinate is 0 or 1)
//         this->singleToCoords(n,local_indices);

//         /* For each dimension, if the local index is 0, use the lower
//            bound of the global index. If the local index is 1, use the
//            "upper" global index.

//            In 2-D:
//            The lower left corner of the element will have local_indices = [0,0];
//            Then, global_indices = [ indices[0], indices[1] ];
//            The upper right corner will have local_indices = [1,1];
//            Then, global_indices = [ indices[0]+1, indices[1]+1 ]; */
//         for( unsigned int d = 0; d < this->m_data.dim(); d++ )
//           {
//             if( local_indices[d] == 0 )
//               global_indices[d] = indices[d];

//             else if( local_indices[d] == 1 )
//               global_indices[d] = indices[d]+1;

//             // This shouldn't happen
//             else
//               queso_error();
//           }

//         /* Now that we have the global indices for each coordinate,
//            we get the "global" index. This is the index into the global
//            values array */
//         unsigned int global = MultiDimensionalIndexing::coordToGlobal( global_indices, this->m_data.get_n_points() );
//         values[n] = this->m_data.get_value(global);
//       }
//   }
  template<class V, class M>
  double GeneralPolynomialChaosInterpolationSurrogate<V,M>::eval_interpolant( const std::vector<double>& x_min,
                                                                      const std::vector<double>& x_max,
                                                                      const std::vector<double>& values,
                                                                      const V & domainVector ) const
  {
    queso_assert_equal_to( x_min.size(), this->m_data.dim() );
    queso_assert_equal_to( x_max.size(), this->m_data.dim() );
    queso_assert_equal_to( values.size(), this->n_coeffs() );
    queso_assert_equal_to( domainVector.sizeGlobal(), this->m_data.dim() );

    double interp_value = 0.0;

    std::vector<unsigned int> indices( this->m_data.dim() );

    for( unsigned int n = 0; n < this->n_coeffs(); n++ )
      {
        this->singleToCoords( n, indices );

        double shape_fn = this->tensor_product_hermite( x_min, x_max, indices, domainVector );

        interp_value += values[n]*shape_fn;
      }

    return interp_value;
  }


  template<class V, class M>
  double GeneralPolynomialChaosInterpolationSurrogate<V,M>::tensor_product_hermite( const std::vector<double>& x_min,
                                                                             const std::vector<double>& x_max,
                                                                             const std::vector<unsigned int>& indices,
                                                                             const V & domainVector ) const
  {
    queso_assert_equal_to( x_min.size(), this->m_data.dim() );
    queso_assert_equal_to( x_max.size(), this->m_data.dim() );
    queso_assert_equal_to( indices.size(), this->m_data.dim() );
    queso_assert_equal_to( domainVector.sizeGlobal(), this->m_data.dim() );
// 
    double value = 1.0;

    for( unsigned int d = 0; d < this->m_data.dim(); d++ )
      {
        value *= this->hermite_poly( x_min[d], x_max[d], domainVector[d], indices[d] );
      }

    return value;
  }

  template<class V, class M>



  template<class V, class M>
  double GeneralPolynomialChaosInterpolationSurrogate<V,M>::hermite_poly( double x0, double x1, double x, unsigned int index , unsigned int order) const
  {
    double value = 0.0;
  	// if (order == 0)
  	
  	// value = (x0+x1/1) * x

   //  for (unsigned int ord = 0; ord < order; order ++)
   //  {

   //  }

   //  // Make sure we're trying to interpolate between the given points
   //  queso_assert_less_equal( x, x1 );
   //  queso_assert_greater_equal( x, x0 );

   //  // Make sure index is either 0 or 1
   //  queso_assert( (index == 0) || (index == 1) );

   //  double value = 0.0;

   //  if( index == 0 )
   //    value = (x-x1)/(x0-x1);
   //  else
   //    value = (x-x0)/(x1-x0);

   //  return value;
  }

  template<class V, class M>
  double GeneralPolynomialChaosInterpolationSurrogate<V,M>::hermite_poly_sym(double x, unsigned int order) const
  {
  	//constructing the symbolic form of Hermite Polynomial
  	//here we differentiate the saying between "order" and "degree". 
  	//"order" is the order of Hermite Polynomial that we are currently calculating
  	//"degree" is the total degree of Hermite polynomial that we are constructing the surrogate for.
  	
  }
  
  template<class V, class M>
  double GeneralPolynomialChaosInterpolationSurrogate<V,M>::hermite_poly_weight(double x0) const
  {
  	//weight function of Hermite Polynomial
  	const double PI  =3.141592653589793238463;
  	double value = (1/std::sqrt(2*PI))*(std::exp(-1*(x0^2)/2));
  	return value;
  }
  double GeneralPolynomialChaosInterpolationSurrogate<V,M>::hermite_poly_weight_gaussian(double x0, int degree) const
  {
  	//weight function of Hermite Polynomial
  	const double PI  =3.141592653589793238463;
  	
  	double value = (1/std::sqrt(2*PI))*(std::exp(-1*(x0^2)/2));
  	return value;
  }

  template<class V, class M>
  double GeneralPolynomialChaosInterpolationSurrogate<V,M>::factorial(unsigned int order) const
  {
  	//factorial for Hermite Polynomial
  	int value = 1;
  	for (int c = 1; c < order +1; c++)
  	{
  		value *= c;
  	}
  	return value;
  }


  template<class V, class M>
  double GeneralPolynomialChaosInterpolationSurrogate<V,M>::hermite_poly_coef_num_int(double x0, double x1, unsigned int order) const
  {
  	//calculating the numerical integration for the Hermite Polynomials coefficient for current order
  	//We differentiate the saying between "order" and "degree". 
  	//"order" is the order of Hermite Polynomial that we are currently calculating
  	//"degree" is the maximum degree of Hermite polynomial that we are constructing the surrogate for.
  	//trapezoidal rule
  	double value = 1;
  	value = (1/2)*((this->hermite_poly_sym(x0, order))*(value(x0)*this->hermite_poly_weight(x0))+(value(x1)this->hermite_poly_sym(x1,order))*(this->hermite_poly_weight(x1)))
  	


  }

  template<class V, class M>
  double GeneralPolynomialChaosInterpolationSurrogate<V,M>::hermite_poly_coef(double x0, double x1, double space, unsigned int order) const
  {
  	//sum over all the numerical integration on current dimentsion to calculate the coefficient of Hermite Polynomial at this dimension of this order
  	double coef = 0;
  	for (double c = x0; c < x1; c+= space)
  	{
  		coef += hermite_poly_coef_num_int(x0, x0 + space, order);

  	}
  	return coef;
  }


} // end namespace QUESO

// Instantiate
template class QUESO::GeneralPolynomialChaosInterpolationSurrogate<QUESO::GslVector,QUESO::GslMatrix>;
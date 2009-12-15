/*--------------------------------------------------------------------------
 *--------------------------------------------------------------------------
 *
 * Copyright (C) 2008 The PECOS Development Team
 *
 * Please see http://pecos.ices.utexas.edu for more information.
 *
 * This file is part of the QUESO Library (Quantification of Uncertainty
 * for Estimation, Simulation and Optimization).
 *
 * QUESO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QUESO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QUESO. If not, see <http://www.gnu.org/licenses/>.
 *
 *--------------------------------------------------------------------------
 *
 * $Id$
 *
 * Brief description of this file: 
 * 
 *--------------------------------------------------------------------------
 *-------------------------------------------------------------------------- */

#ifndef __UQ_GSL_VECTOR_H__
#define __UQ_GSL_VECTOR_H__

#include <uqVector.h>
#include <gsl/gsl_vector.h>

class uqGslVectorClass : public uqVectorClass
{
public:
  uqGslVectorClass();
  uqGslVectorClass(const uqBaseEnvironmentClass& env, const Epetra_Map& map);
  uqGslVectorClass(const uqBaseEnvironmentClass& env, const Epetra_Map& map, double value);
  uqGslVectorClass(const uqBaseEnvironmentClass& env, double d1, double d2, const Epetra_Map& map); // MATLAB linspace
  uqGslVectorClass(const uqGslVectorClass&         v, double d1, double d2);                        // MATLAB linspace
  uqGslVectorClass(const uqGslVectorClass&         y);
 ~uqGslVectorClass();

  uqGslVectorClass& operator= (const uqGslVectorClass& rhs);
  uqGslVectorClass& operator*=(double a);
  uqGslVectorClass& operator/=(double a);
  uqGslVectorClass& operator*=(const uqGslVectorClass& rhs);
  uqGslVectorClass& operator/=(const uqGslVectorClass& rhs);
  uqGslVectorClass& operator+=(const uqGslVectorClass& rhs);
  uqGslVectorClass& operator-=(const uqGslVectorClass& rhs);
            double& operator[](unsigned int i);
      const double& operator[](unsigned int i) const;

  unsigned int sizeLocal        () const;
  unsigned int sizeGlobal       () const;
  double       norm2Sq          () const;
  double       norm2            () const;
  double       norm1            () const;
  double       normInf          () const;
  double       sumOfComponents  () const;
  void         cwSet            (double value);
  void         cwSetGaussian    (const gsl_rng* rng, double mean, double stdDev);
  void         cwSetGaussian    (const gsl_rng* rng, const uqGslVectorClass& meanVec, const uqGslVectorClass& stdDevVec);
  void         cwSetUniform     (const gsl_rng* rng, const uqGslVectorClass& aVec,    const uqGslVectorClass& bVec     );
  void         cwSetInverseGamma(const gsl_rng* rng, const uqGslVectorClass& alpha,   const uqGslVectorClass& beta     );
  void         cwSetConcatenated(const uqGslVectorClass& v1, const uqGslVectorClass& v2);
  void         cwExtract        (unsigned int initialPos, uqGslVectorClass& vec) const;
  void         cwInvert         ();
  void         matlabDiff       (unsigned int firstPositionToStoreDiff, double valueForRemainderPosition, uqGslVectorClass& outputVec) const;
  void         sort             ();
  void         mpiBcast         (int srcRank, const Epetra_MpiComm& bcastComm);
  void         print            (std::ostream& os) const;
  void         subWriteContents (const std::string&            varNamePrefix,
                                 const std::string&            fileName,
                                 const std::set<unsigned int>& allowedSubEnvIds) const;
  void         subReadContents  (const std::string&            fileName,
                                 const std::set<unsigned int>& allowedSubEnvIds);

  bool         atLeastOneComponentSmallerThan(const uqGslVectorClass& rhs) const;
  bool         atLeastOneComponentBiggerThan (const uqGslVectorClass& rhs) const;

  // Necessary for uqGslMatrixClass::invertMultiply() and uqGslMatrixClass::setRow/Column
  gsl_vector*  data                          () const; 

  double       getMaxValue      () const;
  double       getMinValue      () const;
  int          getMaxValueIndex () const;
  int          getMinValueIndex () const;
  void         getMaxValueAndIndex( double& value, int& index );
  void         getMinValueAndIndex( double& value, int& index );

  uqGslVectorClass abs() const;

private:

  void         copy             (const uqGslVectorClass& src);

  gsl_vector* m_vec;
};

uqGslVectorClass operator/    (      double a,            const uqGslVectorClass& x  );
uqGslVectorClass operator/    (const uqGslVectorClass& x, const uqGslVectorClass& y  );
uqGslVectorClass operator*    (      double a,            const uqGslVectorClass& x  );
uqGslVectorClass operator*    (const uqGslVectorClass& x, const uqGslVectorClass& y  );
double           scalarProduct(const uqGslVectorClass& x, const uqGslVectorClass& y  );
uqGslVectorClass operator+    (const uqGslVectorClass& x, const uqGslVectorClass& y  );
uqGslVectorClass operator-    (const uqGslVectorClass& x, const uqGslVectorClass& y  );
std::ostream&    operator<<   (std::ostream& os,          const uqGslVectorClass& obj);

#endif // __UQ_GSL_VECTOR_H__

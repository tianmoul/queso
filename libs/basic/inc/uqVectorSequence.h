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

#ifndef __UQ_VECTOR_SEQUENCE_H__
#define __UQ_VECTOR_SEQUENCE_H__

#include <uqVectorSubset.h>
#include <uqScalarSequence.h>
#include <uqChainStatisticalOptions.h>
#include <uqArrayOfOneDGrids.h>
#include <uqArrayOfOneDTables.h>
#include <uq2dArrayOfStuff.h>
#include <sys/time.h>
#include <fstream>

template <class V, class M>
class uqBaseVectorSequenceClass
{
public:
           uqBaseVectorSequenceClass(const uqVectorSpaceClass<V,M>& vectorSpace,
                                     unsigned int                   sequenceSize,
                                     const std::string&             name);
  virtual ~uqBaseVectorSequenceClass();

  virtual  unsigned int             sequenceSize              () const = 0;
           unsigned int             vectorSize                () const;
  const    uqVectorSpaceClass<V,M>& vectorSpace               () const;
  const    std::string&             name                      () const;
           void                     setName                   (const std::string& newName);
           void                     clear                     ();
  const    V&                       minValues                 () const;
  const    V&                       maxValues                 () const;
  const    V&                       meanValues                () const;
  const    V&                       sampleVarianceValues      () const;
  const    uqBoxSubsetClass<V,M>&   valuesBox                 () const;
           void                     deleteStoredVectors       ();

  virtual  void                     resizeSequence            (unsigned int newSequenceSize) = 0;
  virtual  void                     resetValues               (unsigned int initialPos, unsigned int numPos) = 0;
  virtual  void                     erasePositions            (unsigned int initialPos, unsigned int numPos) = 0;
  virtual  void                     getPositionValues         (unsigned int posId,       V& vec) const = 0;
  virtual  void                     setPositionValues         (unsigned int posId, const V& vec) = 0;
  virtual  void                     setGaussian               (const gsl_rng* rng, const V& meanVec, const V& stdDevVec) = 0;
  virtual  void                     setUniform                (const gsl_rng* rng, const V& aVec,    const V& bVec     ) = 0;
  virtual  void                     uniformlySampledMdf       (const V&                       numEvaluationPointsVec,
                                                               uqArrayOfOneDGridsClass <V,M>& mdfGrids,
                                                               uqArrayOfOneDTablesClass<V,M>& mdfValues) const = 0;
  virtual  void                     uniformlySampledCdf       (const V&                       numEvaluationPointsVec,
                                                               uqArrayOfOneDGridsClass <V,M>& cdfGrids,
                                                               uqArrayOfOneDTablesClass<V,M>& cdfValues) const = 0;
  virtual  void                     unifiedUniformlySampledCdf(const V&                       numEvaluationPointsVec,
                                                               uqArrayOfOneDGridsClass <V,M>& unifiedCdfGrids,
                                                               uqArrayOfOneDTablesClass<V,M>& unifieddfValues) const = 0;

  virtual  void                     mean                      (unsigned int                          initialPos,
                                                               unsigned int                          numPos,
                                                               V&                                    meanVec) const = 0;
  virtual  void                     unifiedMean               (unsigned int                          initialPos,
                                                               unsigned int                          numPos,
                                                               V&                                    unifiedMeanVec) const = 0;
  virtual  void                     sampleVariance            (unsigned int                          initialPos,
                                                               unsigned int                          numPos,
                                                               const V&                              meanVec,
                                                               V&                                    samVec) const = 0;
  virtual  void                     unifiedSampleVariance     (unsigned int                          initialPos,
                                                               unsigned int                          numPos,
                                                               const V&                              unifiedMeanVec,
                                                               V&                                    unifiedSamVec) const = 0;
  virtual  void                     populationVariance        (unsigned int                          initialPos,
                                                               unsigned int                          numPos,
                                                               const V&                              meanVec,
                                                               V&                                    popVec) const = 0;
  virtual  void                     autoCovariance            (unsigned int                          initialPos,
                                                               unsigned int                          numPos,
                                                               const V&                              meanVec,
                                                               unsigned int                          lag,
                                                               V&                                    covVec) const = 0;

  virtual  void                     autoCorrViaDef            (unsigned int                          initialPos,
                                                               unsigned int                          numPos,
                                                               unsigned int                          lag,
                                                               V&                                    corrVec) const = 0;
  virtual  void                     autoCorrViaFft            (unsigned int                          initialPos,
                                                               unsigned int                          numPos,
                                                               const std::vector<unsigned int>&      lags,
                                                               std::vector<V*>&                      corrVecs) const = 0;
  virtual  void                     autoCorrViaFft            (unsigned int                          initialPos,
                                                               unsigned int                          numPos,
                                                               unsigned int                          numSum,
                                                               V&                                    autoCorrsSumVec) const = 0;
  virtual  void                     bmm                       (unsigned int                          initialPos,
                                                               unsigned int                          batchLength,
                                                               V&                                    bmmVec) const = 0;
  virtual  void                     fftForward                (unsigned int                          initialPos,
                                                               unsigned int                          fftSize,
                                                               unsigned int                          paramId,
                                                               std::vector<std::complex<double> >&   fftResult) const = 0;
//virtual  void                     fftInverse                (unsigned int fftSize) = 0;
  virtual  void                     psd                       (unsigned int                          initialPos,
                                                               unsigned int                          numBlocks,
                                                               double                                hopSizeRatio,
                                                               unsigned int                          paramId,
                                                               std::vector<double>&                  psdResult) const = 0;
  virtual  void                     psdAtZero                 (unsigned int                          initialPos,
                                                               unsigned int                          numBlocks,
                                                               double                                hopSizeRatio,
                                                               V&                                    psdVec) const = 0;
  virtual  void                     geweke                    (unsigned int                          initialPos,
                                                               double                                ratioNa,
                                                               double                                ratioNb,
                                                               V&                                    gewVec) const = 0;
  virtual  void                     minMax                    (unsigned int                          initialPos,
                                                               V&                                    minVec,
                                                               V&                                    maxVec) const = 0;
  virtual  void                     unifiedMinMax             (unsigned int                          initialPos,
                                                               V&                                    unifiedMinVec,
                                                               V&                                    unifiedMaxVec) const = 0;
  virtual  void                     histogram                 (unsigned int                          initialPos,
                                                               const V&                              minVec,
                                                               const V&                              maxVec,
                                                               std::vector<V*>&                      centersForAllBins,
                                                               std::vector<V*>&                      quanttsForAllBins) const = 0;
  virtual  void                     unifiedHistogram          (unsigned int                          initialPos,
                                                               const V&                              unifiedMinVec,
                                                               const V&                              unifiedMaxVec,
                                                               std::vector<V*>&                      unifiedCentersForAllBins,
                                                               std::vector<V*>&                      unifiedQuanttsForAllBins) const = 0;
  virtual  void                     interQuantileRange        (unsigned int                          initialPos,
                                                               V&                                    iqrVec) const = 0;
  virtual  void                     unifiedInterQuantileRange (unsigned int                          initialPos,
                                                               V&                                    unifiedIqrVec) const = 0;
  virtual  void                     scalesForKDE              (unsigned int                          initialPos,
                                                               const V&                              iqrVec,
                                                               V&                                    scaleVec) const = 0;
  virtual  void                     unifiedScalesForKDE       (unsigned int                          initialPos,
                                                               const V&                              unifiedIqrVec,
                                                               V&                                    unifiedScaleVec) const = 0;
  virtual  void                     gaussianKDE               (const V&                              evaluationParamVec,
                                                               V&                                    densityVec) const = 0;
  virtual  void                     gaussianKDE               (unsigned int                          initialPos,
                                                               const V&                              scaleVec,
                                                               const std::vector<V*>&                evaluationParamVecs,
                                                               std::vector<V*>&                      densityVecs) const = 0;
  virtual  void                     unifiedGaussianKDE        (unsigned int                          initialPos,
                                                               const V&                              unifiedScaleVec,
                                                               const std::vector<V*>&                unifiedEvaluationParamVecs,
                                                               std::vector<V*>&                      unifiedDensityVecs) const = 0;
  virtual  void                     printContents             (std::ofstream&                        ofsvar) const = 0;
  virtual  void                     printUnifiedContents      (std::ofstream&                        ofsvar) const = 0;
  virtual  void                     printUnifiedContents      (const std::string&                    fileName) const = 0;
  virtual  void                     select                    (const std::vector<unsigned int>&      idsOfUniquePositions) = 0;
  virtual  void                     filter                    (unsigned int                          initialPos,
                                                               unsigned int                          spacing) = 0;

           void                     computeStatistics         (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                               std::ofstream*                        passedOfs);

           void                     computeFilterParams       (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                               std::ofstream*                        passedOfs,
                                                               unsigned int&                         initialPos,
                                                               unsigned int&                         spacing);
protected:
           void                     computeMeanVars           (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                               std::ofstream*                        passedOfs,
                                                               V*                                    meanPtr,
                                                               V*                                    sampleVarPtr,
                                                               V*                                    populVarPtr);
           void                     computeBMM                (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                               const std::vector<unsigned int>&      initialPosForStatistics,
                                                               std::ofstream*                        passedOfs);
           void                     computeFFT                (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                               const std::vector<unsigned int>&      initialPosForStatistics,
                                                               std::ofstream*                        passedOfs);
           void                     computePSD                (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                               const std::vector<unsigned int>&      initialPosForStatistics,
                                                               std::ofstream*                        passedOfs);
           void                     computePSDAtZero          (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                               const std::vector<unsigned int>&      initialPosForStatistics,
                                                               std::ofstream*                        passedOfs);
           void                     computeGeweke             (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                               const std::vector<unsigned int>&      initialPosForStatistics,
                                                               std::ofstream*                        passedOfs);
           void                     computeAutoCorrViaDef     (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                               const std::vector<unsigned int>&      initialPosForStatistics,
                                                               const std::vector<unsigned int>&      lagsForCorrs,
                                                               std::ofstream*                        passedOfs);
           void                     computeAutoCorrViaFFT     (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                               const std::vector<unsigned int>&      initialPosForStatistics,
                                                               const std::vector<unsigned int>&      lagsForCorrs,
                                                               std::ofstream*                        passedOfs);
           void                     computeHistKde            (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                               std::ofstream*                        passedOfs);
           void                     computeCovCorrMatrices    (const uqChainStatisticalOptionsClass& statisticalOptions,
                                                               std::ofstream*                        passedOfs);

  virtual  void                     extractScalarSeq          (unsigned int                          initialPos,
                                                               unsigned int                          spacing,
                                                               unsigned int                          numPos,
                                                               unsigned int                          paramId,
                                                               uqScalarSequenceClass<double>&        scalarSeq) const = 0;
  virtual  void                     extractRawData            (unsigned int                          initialPos,
                                                               unsigned int                          spacing,
                                                               unsigned int                          numPos,
                                                               unsigned int                          paramId,
                                                               std::vector<double>&                  rawData) const = 0;

  const uqBaseEnvironmentClass&  m_env;
  const uqVectorSpaceClass<V,M>& m_vectorSpace;
  std::string                    m_name;

  mutable uqFftClass<double>*    m_fftObj;
  mutable V*                     m_minValues;
  mutable V*                     m_maxValues;
  mutable V*                     m_meanValues;
  mutable V*                     m_sampleVarianceValues;
  mutable uqBoxSubsetClass<V,M>* m_valuesBox;
};

template <class V, class M>
uqBaseVectorSequenceClass<V,M>::uqBaseVectorSequenceClass(
  const uqVectorSpaceClass<V,M>& vectorSpace,
  unsigned int                   sequenceSize,
  const std::string&             name)
  :
  m_env                 (vectorSpace.env()),
  m_vectorSpace         (vectorSpace),
  m_name                (name),
  m_fftObj              (new uqFftClass<double>(m_env)),
  m_minValues           (NULL),
  m_maxValues           (NULL),
  m_meanValues          (NULL),
  m_sampleVarianceValues(NULL),
  m_valuesBox           (NULL)
{
}

template <class V, class M>
uqBaseVectorSequenceClass<V,M>::~uqBaseVectorSequenceClass()
{
  //clear();
  if (m_valuesBox           ) delete m_valuesBox;
  if (m_sampleVarianceValues) delete m_sampleVarianceValues;
  if (m_meanValues          ) delete m_meanValues;
  if (m_maxValues           ) delete m_maxValues;
  if (m_minValues           ) delete m_minValues;
  if (m_fftObj != NULL      ) delete m_fftObj;
}

template <class V, class M>
unsigned int
uqBaseVectorSequenceClass<V,M>::vectorSize() const
{
  return m_vectorSpace.dim();
}

template <class V, class M>
const uqVectorSpaceClass<V,M>&
uqBaseVectorSequenceClass<V,M>::vectorSpace() const
{
  return m_vectorSpace;
}

template <class V, class M>
const std::string&
uqBaseVectorSequenceClass<V,M>::name() const
{
  return m_name;
}

template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::setName(const std::string& newName)
{
  m_name = newName;
  return;
}

template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::clear()
{
  unsigned int numPos = this->sequenceSize();
  if (numPos) {
    this->resetValues(0,numPos);
    this->resizeSequence(0);
  }

 return;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::minValues() const
{
  if (m_minValues == NULL) {
    m_minValues = m_vectorSpace.newVector();
    if (m_maxValues == NULL) m_maxValues = m_vectorSpace.newVector();
    minMax(0,*m_minValues,*m_maxValues);
  }

  return *m_minValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::maxValues() const
{
  if (m_maxValues == NULL) {
    if (m_minValues == NULL) m_minValues = m_vectorSpace.newVector();
    m_maxValues = m_vectorSpace.newVector();
    minMax(0,*m_minValues,*m_maxValues);
  }

  return *m_maxValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::meanValues() const
{
  if (m_meanValues == NULL) {
    m_meanValues = m_vectorSpace.newVector();
    mean(0,sequenceSize(),*m_meanValues);
  }

  return *m_meanValues;
}

template <class V, class M>
const V&
uqBaseVectorSequenceClass<V,M>::sampleVarianceValues() const
{
  if (m_sampleVarianceValues == NULL) {
    m_sampleVarianceValues = m_vectorSpace.newVector();
    sampleVariance(0,sequenceSize(),meanValues(),*m_sampleVarianceValues);
  }

  return *m_sampleVarianceValues;
}

template <class V, class M>
const uqBoxSubsetClass<V,M>&
uqBaseVectorSequenceClass<V,M>::valuesBox() const
{
  if (m_valuesBox == NULL) {
    m_valuesBox = new uqBoxSubsetClass<V,M>(m_name.c_str(),
                                            m_vectorSpace,
                                            this->minValues(),
                                            this->maxValues());
  }

  return *m_valuesBox;
}

template <class V, class M>
void
uqBaseVectorSequenceClass<V,M>::deleteStoredVectors()
{
  delete m_minValues;
  delete m_maxValues;
  delete m_meanValues;
  delete m_sampleVarianceValues;

  m_minValues            = NULL;
  m_maxValues            = NULL;
  m_meanValues           = NULL;
  m_sampleVarianceValues = NULL;

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeStatistics(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                        passedOfs)
{
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n"
                           << "\n-----------------------------------------------------"
                           << "\n Computing statistics for chain " << m_name << " ..."
                           << "\n-----------------------------------------------------"
                           << "\n"
                           << std::endl;
  }

  bool okSituation = ((passedOfs == NULL                          ) ||
                      (passedOfs != NULL) && (m_env.subRank() >= 0));
  UQ_FATAL_TEST_MACRO(!okSituation,
                      m_env.rank(),
                      "uqBaseVectorSequenceClass<V,M>::computeStatistics()",
                      "unexpected combination of file pointer and subRank");

  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  // Set initial positions for the computation of chain statistics
  std::vector<unsigned int> initialPosForStatistics(statisticalOptions.initialDiscardedPortions().size(),0);
  for (unsigned int i = 0; i < initialPosForStatistics.size(); ++i) {
    initialPosForStatistics[i] = (unsigned int) (statisticalOptions.initialDiscardedPortions()[i] * (double) this->sequenceSize());
  }
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "In uqBaseVectorSequenceClass<V,M>::computeStatistics(): initial positions for statistics =";
    for (unsigned int i = 0; i < initialPosForStatistics.size(); ++i) {
      *m_env.subScreenFile() << " " << initialPosForStatistics[i];
    }
    *m_env.subScreenFile() << std::endl;
  }

  //****************************************************
  // Compute mean, sample std, population std
  //****************************************************
  this->computeMeanVars(statisticalOptions,
                        passedOfs,
                        NULL,
                        NULL,
                        NULL);

  //****************************************************
  // Compute variance of sample mean through the 'batch means method' (BMM)
  //****************************************************
  if ((statisticalOptions.bmmRun()               ) &&
      (initialPosForStatistics.size()         > 0) &&
      (statisticalOptions.bmmLengths().size() > 0)) { 
    this->computeBMM(statisticalOptions,
                     initialPosForStatistics,
                     passedOfs);
  }

  //****************************************************
  // Compute FFT of chain, for one parameter only
  //****************************************************
  if ((statisticalOptions.fftCompute()   ) &&
      (initialPosForStatistics.size() > 0)) {
    this->computeFFT(statisticalOptions,
                     initialPosForStatistics,
                     passedOfs);
  }

  //****************************************************
  // Compute power spectral density (PSD) of chain, for one parameter only
  //****************************************************
  if ((statisticalOptions.psdCompute()   ) &&
      (initialPosForStatistics.size() > 0)) {
    this->computePSD(statisticalOptions,
                     initialPosForStatistics,
                     passedOfs);
  }

  //****************************************************
  // Compute power spectral density (PSD) of chain at zero frequency
  //****************************************************
  if ((statisticalOptions.psdAtZeroCompute()             ) &&
      (initialPosForStatistics.size()                 > 0) &&
      (statisticalOptions.psdAtZeroNumBlocks().size() > 0)) { 
    this->computePSDAtZero(statisticalOptions,
                           initialPosForStatistics,
                           passedOfs);
  }

  //****************************************************
  // Compute Geweke
  //****************************************************
  if ((statisticalOptions.gewekeCompute()) &&
      (initialPosForStatistics.size() > 0)) {
    this->computeGeweke(statisticalOptions,
                        initialPosForStatistics,
                        passedOfs);
  }

  // Set lags for the computation of chain autocorrelations
  std::vector<unsigned int> lagsForCorrs(statisticalOptions.autoCorrNumLags(),1);
  for (unsigned int i = 1; i < lagsForCorrs.size(); ++i) {
    lagsForCorrs[i] = statisticalOptions.autoCorrSecondLag() + (i-1)*statisticalOptions.autoCorrLagSpacing();
  }

  //****************************************************
  // Compute autocorrelation coefficients via definition
  //****************************************************
  if ((statisticalOptions.autoCorrComputeViaDef()) &&
      (initialPosForStatistics.size() > 0    ) &&
      (lagsForCorrs.size()            > 0    )) { 
    this->computeAutoCorrViaDef(statisticalOptions,
                                initialPosForStatistics,
                                lagsForCorrs,
                                passedOfs);
  }

  //****************************************************
  // Compute autocorrelation coefficients via FFT
  //****************************************************
  if ((statisticalOptions.autoCorrComputeViaFft()) &&
      (initialPosForStatistics.size() > 0    ) &&
      (lagsForCorrs.size()            > 0    )) { 
    this->computeAutoCorrViaFFT(statisticalOptions,
                                initialPosForStatistics,
                                lagsForCorrs,
                                passedOfs);
  }

  //****************************************************
  // Compute histogram and/or Kde
  //****************************************************
  if ((statisticalOptions.histCompute()) ||
      (statisticalOptions.kdeCompute() )) {
    this->computeHistKde(statisticalOptions,
                         passedOfs);
  }

  //****************************************************
  // Compute covariance and correlation matrices
  //****************************************************
  if ((statisticalOptions.covMatrixCompute ()) ||
      (statisticalOptions.corrMatrixCompute())) {
    this->computeCovCorrMatrices(statisticalOptions,
                                 passedOfs);
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "All statistics took " << tmpRunTime
                           << " seconds"
                           << std::endl;
  }

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\n Finished computing statistics for chain " << m_name
                           << "\n-----------------------------------------------------"
                           << "\n"
                           << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeMeanVars(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                        passedOfs,
  V*                                    meanPtr,
  V*                                    sampleVarPtr,
  V*                                    populVarPtr)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\nComputing mean, sample variance and population variance"
                           << std::endl;
  }

  V chainMean(m_vectorSpace.zeroVector());
  this->mean(0,
             this->sequenceSize(),
             chainMean);

  V chainSampleVariance(m_vectorSpace.zeroVector());
  this->sampleVariance(0,
                       this->sequenceSize(),
                       chainMean,
                       chainSampleVariance);

  if ((m_env.verbosity() >= 5) && (m_env.subScreenFile())) {
    *m_env.subScreenFile() << "In uqBaseVectorSequenceClass<V,M>::computeMeanVars()"
                           << ": chainMean.size() = "           << chainMean.size()
                           << ", chainMean = "                  << chainMean
                           << ", chainSampleVariance.size() = " << chainSampleVariance.size()
                           << ", chainSampleVariance = "        << chainSampleVariance
                           << std::endl;
  }

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\nEstimated variance of sample mean for the whole chain " << m_name
                           << ", under independence assumption:"
                           << std::endl;
  }
  V estimatedVarianceOfSampleMean(chainSampleVariance);
  estimatedVarianceOfSampleMean /= (double) this->sequenceSize();
  bool savedVectorPrintState = estimatedVarianceOfSampleMean.getPrintHorizontally();
  estimatedVarianceOfSampleMean.setPrintHorizontally(false);
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << estimatedVarianceOfSampleMean
                           << std::endl;
  }
  estimatedVarianceOfSampleMean.setPrintHorizontally(savedVectorPrintState);

  V chainPopulationVariance(m_vectorSpace.zeroVector());
  this->populationVariance(0,
                           this->sequenceSize(),
                           chainMean,
                           chainPopulationVariance);

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Mean and variances took " << tmpRunTime
                           << " seconds"
                           << std::endl;
  }

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\nMean, sample std, population std"
                           << std::endl;
    char line[512];
    sprintf(line,"%s%4s%s%9s%s%9s%s",
	    "Parameter",
            " ",
            "Mean",
            " ",
            "SampleStd",
            " ",
            "Popul.Std");
    *m_env.subScreenFile() << line;

    for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
      sprintf(line,"\n%8.8s%2s%11.4e%2s%11.4e%2s%11.4e",
              m_vectorSpace.componentName(i).c_str(), /*.*/
              " ",
	      chainMean[i],
              " ",
              sqrt(chainSampleVariance[i]),
              " ",
              sqrt(chainPopulationVariance[i]));
      *m_env.subScreenFile() << line;
    }
    *m_env.subScreenFile() << std::endl;
  }

  if (meanPtr     ) *meanPtr      = chainMean;
  if (sampleVarPtr) *sampleVarPtr = chainSampleVariance;
  if (populVarPtr ) *populVarPtr  = chainPopulationVariance;

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeBMM(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\nComputing variance of sample mean through BMM"
                           << std::endl;
  }

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "In uqBaseVectorSequenceClass<V,M>::computeBMM(): lengths for batchs in BMM =";
    for (unsigned int i = 0; i < statisticalOptions.bmmLengths().size(); ++i) {
      *m_env.subScreenFile() << " " << statisticalOptions.bmmLengths()[i];
    }
    *m_env.subScreenFile() << std::endl;
  }

  uq2dArrayOfStuff<V> _2dArrayOfBMM(initialPosForStatistics.size(),statisticalOptions.bmmLengths().size());
  for (unsigned int i = 0; i < _2dArrayOfBMM.numRows(); ++i) {
    for (unsigned int j = 0; j < _2dArrayOfBMM.numCols(); ++j) {
      _2dArrayOfBMM.setLocation(i,j,new V(m_vectorSpace.zeroVector()) /*.*/);
    }
  }
  V bmmVec(m_vectorSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPos = initialPosForStatistics[initialPosId];
    for (unsigned int batchLengthId = 0; batchLengthId < statisticalOptions.bmmLengths().size(); batchLengthId++) {
      unsigned int batchLength = statisticalOptions.bmmLengths()[batchLengthId];
      this->bmm(initialPos,
                batchLength,
                bmmVec);
      _2dArrayOfBMM(initialPosId,batchLengthId) = bmmVec;
    }
  }

  if (m_env.subScreenFile()) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      *m_env.subScreenFile() << "\nEstimated variance of sample mean, through batch means method, for subchain beggining at position " << initialPosForStatistics[initialPosId]
                             << " (each column corresponds to a batch length)"
                             << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      *m_env.subScreenFile() << line;
      for (unsigned int batchLengthId = 0; batchLengthId < statisticalOptions.bmmLengths().size(); batchLengthId++) {
        sprintf(line,"%10s%3d",
                " ",
                statisticalOptions.bmmLengths()[batchLengthId]);
        *m_env.subScreenFile() << line;
      }

      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.componentName(i).c_str() /*.*/);
        *m_env.subScreenFile() << line;
        for (unsigned int batchLengthId = 0; batchLengthId < statisticalOptions.bmmLengths().size(); batchLengthId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfBMM(initialPosId,batchLengthId)[i]);
          *m_env.subScreenFile() << line;
        }
      }
      *m_env.subScreenFile() << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Chain BMM took " << tmpRunTime
                           << " seconds"
                           << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeFFT(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\nComputing FFT of chain on parameter of id = " << statisticalOptions.fftParamId()
                           << std::endl;
  }

  std::vector<std::complex<double> > forwardResult(0,std::complex<double>(0.,0.));
  std::vector<std::complex<double> > inverseResult(0,std::complex<double>(0.,0.));
  uqFftClass<std::complex<double> > fftObj(m_env);
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    this->fftForward(initialPosition,
                     statisticalOptions.fftSize(),
                     statisticalOptions.fftParamId(),
                     forwardResult);

    if (statisticalOptions.fftWrite() && passedOfs) {
      std::ofstream& ofsvar = *passedOfs;
      ofsvar << m_name << "_fftInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << " = zeros(" << 1
             << ","                                                                                                              << forwardResult.size()
             << ");"
             << std::endl;
      for (unsigned int j = 0; j < forwardResult.size(); ++j) {
        ofsvar << m_name << "_fftInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << "(" << 1
               << ","                                                                                               << j+1
               << ") = "                                                                                            << forwardResult[j].real()
               << " + i*"                                                                                           << forwardResult[j].imag()
               << ";"
               << std::endl;
      }
    } // if write

    if (statisticalOptions.fftTestInversion()) {
      fftObj.inverse(forwardResult,
                     statisticalOptions.fftSize(),
                     inverseResult);
      if (statisticalOptions.fftWrite() && passedOfs) {
        std::ofstream& ofsvar = *passedOfs;
        ofsvar << m_name << "_iftInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << " = zeros(" << 1
               << ","                                                                                                              << inverseResult.size()
               << ");"
               << std::endl;
        for (unsigned int j = 0; j < inverseResult.size(); ++j) {
          ofsvar << m_name << "_iftInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << "(" << 1
                 << ","                                                                                               << j+1
                 << ") = "                                                                                            << inverseResult[j].real()
                 << " + i*"                                                                                           << inverseResult[j].imag()
                 << ";"
                 << std::endl;
        }
      } // if write
    }
  } // for initialPosId

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Chain FFT took " << tmpRunTime
                           << " seconds"
                           << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computePSD(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\nComputing PSD of chain on parameter of id = " << statisticalOptions.psdParamId()
                           << std::endl;
  }

  std::vector<double> psdResult(0,0.);
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    this->psd(initialPosition,
              statisticalOptions.psdNumBlocks(),
              statisticalOptions.psdHopSizeRatio(),
              statisticalOptions.psdParamId(),
              psdResult);

    if (statisticalOptions.psdWrite() && passedOfs) {
      std::ofstream& ofsvar = *passedOfs;
      ofsvar << m_name << "_psdInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << " = zeros(" << 1
             << ","                                                                                                              << psdResult.size()
             << ");"
             << std::endl;
      for (unsigned int j = 0; j < psdResult.size(); ++j) {
        ofsvar << m_name << "_psdInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << "(" << 1
               << ","                                                                                                      << j+1
               << ") = "                                                                                                   << psdResult[j]
               << ";"
               << std::endl;
      }
    } // if write
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Chain PSD took " << tmpRunTime
                           << " seconds"
                           << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computePSDAtZero(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\nComputing PSD at frequency zero for all parameters"
                           << std::endl;
  }

  uq2dArrayOfStuff<V> _2dArrayOfPSDAtZero(initialPosForStatistics.size(),statisticalOptions.psdAtZeroNumBlocks().size());
  for (unsigned int i = 0; i < _2dArrayOfPSDAtZero.numRows(); ++i) {
    for (unsigned int j = 0; j < _2dArrayOfPSDAtZero.numCols(); ++j) {
      _2dArrayOfPSDAtZero.setLocation(i,j,new V(m_vectorSpace.zeroVector()) /*.*/);
    }
  }
  V psdVec(m_vectorSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
      unsigned int numBlocks = statisticalOptions.psdAtZeroNumBlocks()[numBlocksId];
      this->psdAtZero(initialPosition,
                      numBlocks,
                      statisticalOptions.psdAtZeroHopSizeRatio(),
                      psdVec);
      _2dArrayOfPSDAtZero(initialPosId,numBlocksId) = psdVec;
    }
  }

  // Display PSD at frequency zero
  if ((statisticalOptions.psdAtZeroDisplay()) && (m_env.subScreenFile())) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];
      *m_env.subScreenFile() << "\nComputed PSD at frequency zero for subchain beggining at position " << initialPos
                             << ", so effective data size = " << this->sequenceSize() - initialPos
                             << " (each column corresponds to a number of blocks)"
                             << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      *m_env.subScreenFile() << line;
      for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
        sprintf(line,"%10s%3d",
                " ",
                statisticalOptions.psdAtZeroNumBlocks()[numBlocksId]);
        *m_env.subScreenFile() << line;
      }

      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.componentName(i).c_str() /*.*/);
        *m_env.subScreenFile() << line;
        for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]);
          *m_env.subScreenFile() << line;
        }
      }
      *m_env.subScreenFile() << std::endl;
    }
  }

  // Display estimated variance of sample mean through PSD
  if (/*(statisticalOptions.psdAtZeroDisplay()) &&*/ (m_env.subScreenFile())) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];
      *m_env.subScreenFile() << "\nEstimated variance of sample mean, through psd, for subchain beggining at position " << initialPos
                << ", so effective data size = " << this->sequenceSize() - initialPos
                << " (each column corresponds to a number of blocks)"
                << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      *m_env.subScreenFile() << line;
      for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
        sprintf(line,"%10s%3d",
                " ",
                statisticalOptions.psdAtZeroNumBlocks()[numBlocksId]);
        *m_env.subScreenFile() << line;
      }

      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.componentName(i).c_str() /*.*/);
        *m_env.subScreenFile() << line;
        for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  2.*M_PI*_2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]/(double) (this->sequenceSize() - initialPos));
          *m_env.subScreenFile() << line;
        }
      }
      *m_env.subScreenFile() << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Chain PSD at frequency zero took " << tmpRunTime 
                           << " seconds"
                           << std::endl;
  }

  // Write PSD at frequency zero
  if (statisticalOptions.psdAtZeroWrite() && passedOfs) {
    std::ofstream& ofsvar = *passedOfs;
    ofsvar << m_name << "_psdAtZeroNumBlocks_sub" << m_env.subIdString() << " = zeros(" << 1
           << ","                                                                       << statisticalOptions.psdAtZeroNumBlocks().size()
           << ");"
           << std::endl;
    for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
      ofsvar << m_name << "_psdAtZeroNumBlocks_sub" << m_env.subIdString() << "(" << 1
             << ","                                                               << numBlocksId+1
             << ") = "                                                            << statisticalOptions.psdAtZeroNumBlocks()[numBlocksId]
             << ";"
             << std::endl;
    }

    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      ofsvar << m_name << "_psdAtZeroInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << " = zeros(" << this->vectorSize() /*.*/
             << ","                                                                                                                    << statisticalOptions.psdAtZeroNumBlocks().size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int numBlocksId = 0; numBlocksId < statisticalOptions.psdAtZeroNumBlocks().size(); numBlocksId++) {
          ofsvar << m_name << "_psdAtZeroInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << "(" << i+1
                 << ","                                                                                                            << numBlocksId+1
                 << ") = "                                                                                                         << _2dArrayOfPSDAtZero(initialPosId,numBlocksId)[i]
                 << ";"
                 << std::endl;
        }
      }
    }
  } 

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeGeweke(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\nComputing Geweke coefficients"
                           << std::endl;
  }

  std::vector<V*> vectorOfGeweke(initialPosForStatistics.size(),NULL);
  V gewVec(m_vectorSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPosition = initialPosForStatistics[initialPosId];
    this->geweke(initialPosition,
                 statisticalOptions.gewekeNaRatio(),
                 statisticalOptions.gewekeNbRatio(),
                 gewVec);
    vectorOfGeweke[initialPosId] = new V(gewVec);
  }

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\nComputed Geweke coefficients with 10% and 50% percentages"
                           << " (each column corresponds to a different initial position on the full chain)"
                           << std::endl;

    char line[512];
    sprintf(line,"%s",
            "Parameter");
    *m_env.subScreenFile() << line;
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      sprintf(line,"%10s%3d",
              " ",
              initialPosForStatistics[initialPosId]);
      *m_env.subScreenFile() << line;
    }

    for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
      sprintf(line,"\n%9.9s",
              m_vectorSpace.componentName(i).c_str() /*.*/);
      *m_env.subScreenFile() << line;
      for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
        sprintf(line,"%2s%11.4e",
                " ",
                (*(vectorOfGeweke[initialPosId]))[i]);
        *m_env.subScreenFile() << line;
      }
    }
    *m_env.subScreenFile() << std::endl;
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Chain Geweke took " << tmpRunTime
                           << " seconds"
                           << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeAutoCorrViaDef(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::vector<unsigned int>&      lagsForCorrs,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\nComputing autocorrelation coefficients (via def)"
                           << std::endl;
  }

  if (statisticalOptions.autoCorrDisplay() && (m_env.subScreenFile())) {
    *m_env.subScreenFile() << "In uqBaseVectorSequenceClass<V,M>::computeAutoCorrViaDef(): lags for autocorrelation (via def) = ";
    for (unsigned int i = 0; i < lagsForCorrs.size(); ++i) {
      *m_env.subScreenFile() << " " << lagsForCorrs[i];
    }
    *m_env.subScreenFile() << std::endl;
  }

  uq2dArrayOfStuff<V> _2dArrayOfAutoCorrs(initialPosForStatistics.size(),lagsForCorrs.size());
  for (unsigned int i = 0; i < _2dArrayOfAutoCorrs.numRows(); ++i) {
    for (unsigned int j = 0; j < _2dArrayOfAutoCorrs.numCols(); ++j) {
      _2dArrayOfAutoCorrs.setLocation(i,j,new V(m_vectorSpace.zeroVector()) /*.*/);
    }
  }
  //V corrVec(m_vectorSpace.zeroVector());
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    unsigned int initialPos = initialPosForStatistics[initialPosId];
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      unsigned int lag = lagsForCorrs[lagId];
      this->autoCorrViaDef(initialPos,
                           this->sequenceSize()-initialPos,
                           lag,
                           _2dArrayOfAutoCorrs(initialPosId,lagId));
      //_2dArrayOfAutoCorrs(initialPosId,lagId) = corrVec;
    }
  }

  // It is not practical to compute the variance of sample mean by computing the autocorrelations via definition for each lag
  // The code computes the variance of sample mean by computing the autocorrelations via fft, below, in another routine

  if ((statisticalOptions.autoCorrDisplay()) && (m_env.subScreenFile())) {
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      *m_env.subScreenFile() << "\nComputed autocorrelation coefficients (via def), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                             << " (each column corresponds to a different lag)"
                             << std::endl;
      char line[512];
      sprintf(line,"%s",
              "Parameter");
      *m_env.subScreenFile() << line;
      for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
        sprintf(line,"%10s%3d",
                " ",
                lagsForCorrs[lagId]);
        *m_env.subScreenFile() << line;
      }

      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        sprintf(line,"\n%9.9s",
                m_vectorSpace.componentName(i).c_str() /*.*/);
        *m_env.subScreenFile() << line;
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          sprintf(line,"%2s%11.4e",
                  " ",
                  _2dArrayOfAutoCorrs(initialPosId,lagId)[i]);
          *m_env.subScreenFile() << line;
        }
      }
      *m_env.subScreenFile() << std::endl;
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Chain autocorrelation (via def) took " << tmpRunTime
                           << " seconds"
                           << std::endl;
  }

  // Write autocorrelations
  if (statisticalOptions.autoCorrWrite() && passedOfs) {
    std::ofstream& ofsvar = *passedOfs;
    ofsvar << m_name << "_corrViaDefLags_sub" << m_env.subIdString() << " = zeros(" << 1
           << ","                                                                   << lagsForCorrs.size()
           << ");"
           << std::endl;
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      ofsvar << m_name << "_corrViaDefLags_sub" << m_env.subIdString() << "(" << 1
             << ","                                                           << lagId+1
             << ") = "                                                        << lagsForCorrs[lagId]
             << ";"
             << std::endl;
    }

    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      ofsvar << m_name << "_corrViaDefInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << " = zeros(" << this->vectorSize() /*.*/
             << ","                                                                                                                     << lagsForCorrs.size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          ofsvar << m_name << "_corrViaDefInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << "(" << i+1
                 << ","                                                                                                             << lagId+1
                 << ") = "                                                                                                          << _2dArrayOfAutoCorrs(initialPosId,lagId)[i]
                 << ";"
                 << std::endl;
        }
      }
    }
  } 

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeAutoCorrViaFFT(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  const std::vector<unsigned int>&      initialPosForStatistics,
  const std::vector<unsigned int>&      lagsForCorrs,
  std::ofstream*                        passedOfs)
{
  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;
  iRC = gettimeofday(&timevalTmp, NULL);
  double tmpRunTime = 0.;

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\nComputing autocorrelation coefficients (via fft)"
                           << std::endl;
  }

  if (statisticalOptions.autoCorrDisplay() && (m_env.subScreenFile())) {
    *m_env.subScreenFile() << "In uqBaseVectorSequenceClass<V,M>::computeAutoCorrViaFFT(): lags for autocorrelation (via fft) = ";
    for (unsigned int i = 0; i < lagsForCorrs.size(); ++i) {
      *m_env.subScreenFile() << " " << lagsForCorrs[i];
     }
     *m_env.subScreenFile() << std::endl;
  }

  uq2dArrayOfStuff<V> _2dArrayOfAutoCorrs(initialPosForStatistics.size(),lagsForCorrs.size());
  for (unsigned int i = 0; i < _2dArrayOfAutoCorrs.numRows(); ++i) {
    for (unsigned int j = 0; j < _2dArrayOfAutoCorrs.numCols(); ++j) {
      _2dArrayOfAutoCorrs.setLocation(i,j,new V(m_vectorSpace.zeroVector()) /*.*/);
    }
  }
  std::vector<V*> corrVecs(lagsForCorrs.size(),NULL);
  std::vector<V*> corrSumVecs(initialPosForStatistics.size(),NULL);
  for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
    corrSumVecs[initialPosId] = new V(m_vectorSpace.zeroVector()) /*.*/;
    unsigned int initialPos = initialPosForStatistics[initialPosId];
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      corrVecs[lagId] = new V(m_vectorSpace.zeroVector()) /*.*/;
    }
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "In uqBaseVectorSequenceClass<V,M>::computeAutoCorrViaFFT()"
                             << ": about to call chain.autoCorrViaFft()"
                             << " with initialPos = "      << initialPos
                             << ", numPos = "              << this->sequenceSize()-initialPos
                             << ", lagsForCorrs.size() = " << lagsForCorrs.size()
                             << ", corrVecs.size() = "     << corrVecs.size()
                             << std::endl;
    }
    this->autoCorrViaFft(initialPos,
                         this->sequenceSize()-initialPos, // Use all possible data positions
                         lagsForCorrs,
                         corrVecs);
    this->autoCorrViaFft(initialPos,
                         this->sequenceSize()-initialPos, // Use all possible data positions
                         (unsigned int) (1.0 * (double) (this->sequenceSize()-initialPos)), // CHECK
                         *corrSumVecs[initialPosId]); // Sum of all possibly computable autocorrelations, not only the asked ones in lagsForCorrs
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      _2dArrayOfAutoCorrs(initialPosId,lagId) = *(corrVecs[lagId]);
    }
  }
  for (unsigned int j = 0; j < corrVecs.size(); ++j) {
    if (corrVecs[j] != NULL) delete corrVecs[j];
  }

  if (statisticalOptions.autoCorrDisplay()) {
    V chainMean                    (m_vectorSpace.zeroVector());
    V chainSampleVariance          (m_vectorSpace.zeroVector());
    V estimatedVarianceOfSampleMean(m_vectorSpace.zeroVector());
    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      unsigned int initialPos = initialPosForStatistics[initialPosId];

      this->mean(initialPos,
                 this->sequenceSize()-initialPos,
                 chainMean);

      this->sampleVariance(initialPos,
                           this->sequenceSize()-initialPos,
                           chainMean,
                           chainSampleVariance);

      if (m_env.subScreenFile()) {
        *m_env.subScreenFile() << "\nEstimated variance of sample mean, through autocorrelation (via fft), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                               << std::endl;
      }
      estimatedVarianceOfSampleMean.cwSet(-1.); // Yes, '-1' because the autocorrelation at lag 0, which values '+1', is already counted in the sum
      estimatedVarianceOfSampleMean += 2.* (*corrSumVecs[initialPosId]);
      estimatedVarianceOfSampleMean *= chainSampleVariance;
      estimatedVarianceOfSampleMean /= (double) (this->sequenceSize() - initialPos);
      bool savedVectorPrintState = estimatedVarianceOfSampleMean.getPrintHorizontally();
      estimatedVarianceOfSampleMean.setPrintHorizontally(false);
      if (m_env.subScreenFile()) {
        *m_env.subScreenFile() << estimatedVarianceOfSampleMean
                               << std::endl;
      }
      estimatedVarianceOfSampleMean.setPrintHorizontally(savedVectorPrintState);

      if (m_env.subScreenFile()) {
        *m_env.subScreenFile() << "\nComputed autocorrelation coefficients (via fft), for subchain beggining at position " << initialPosForStatistics[initialPosId]
                               << " (each column corresponds to a different lag)"
                               << std::endl;

        char line[512];
        sprintf(line,"%s",
                "Parameter");
        *m_env.subScreenFile() << line;
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          sprintf(line,"%10s%3d",
                  " ",
                  lagsForCorrs[lagId]);
          *m_env.subScreenFile() << line;
        }

        for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
          sprintf(line,"\n%9.9s",
                  m_vectorSpace.componentName(i).c_str() /*.*/);
          *m_env.subScreenFile() << line;
          for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
            sprintf(line,"%2s%11.4e",
                    " ",
                    _2dArrayOfAutoCorrs(initialPosId,lagId)[i]);
            *m_env.subScreenFile() << line;
          }
        }
        *m_env.subScreenFile() << std::endl;
      }
    }
  }

  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Chain autocorrelation (via fft) took " << tmpRunTime
                           << " seconds"
                           << std::endl;
  }

  // Write autocorrelations
  if (statisticalOptions.autoCorrWrite() && passedOfs) {
    std::ofstream& ofsvar = *passedOfs;
    ofsvar << m_name << "_corrViaFftLags_sub" << m_env.subIdString() << " = zeros(" << 1
           << ","                                                                   << lagsForCorrs.size()
           << ");"
           << std::endl;
    for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
      ofsvar << m_name << "_corrViaFftLags_sub" << m_env.subIdString() << "(" << 1
             << ","                                                           << lagId+1
             << ") = "                                                        << lagsForCorrs[lagId]
             << ";"
             << std::endl;
    }

    for (unsigned int initialPosId = 0; initialPosId < initialPosForStatistics.size(); initialPosId++) {
      ofsvar << m_name << "_corrViaFftInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << " = zeros(" << this->vectorSize() /*.*/
             << ","                                                                                                                     << lagsForCorrs.size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int lagId = 0; lagId < lagsForCorrs.size(); lagId++) {
          ofsvar << m_name << "_corrViaFftInitPos" << initialPosForStatistics[initialPosId] << "_sub" << m_env.subIdString() << "(" << i+1
                 << ","                                                                                                             << lagId+1
                 << ") = "                                                                                                          << _2dArrayOfAutoCorrs(initialPosId,lagId)[i]
                 << ";"
                 << std::endl;
        }
      }
    }
  } 

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeFilterParams(
  const uqChainStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                        passedOfs,
  unsigned int&                         initialPos,
  unsigned int&                         spacing)
{
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n"
                           << "\n-----------------------------------------------------"
                           << "\n Computing filter parameters for chain " << m_name << " ..."
                           << "\n-----------------------------------------------------"
                           << "\n"
                           << std::endl;
  }

  bool okSituation = ((passedOfs == NULL                          ) ||
                      (passedOfs != NULL) && (m_env.subRank() >= 0));
  UQ_FATAL_TEST_MACRO(!okSituation,
                      m_env.rank(),
                      "uqBaseVectorSequenceClass<V,M>::computeFilterParams()",
                      "unexpected combination of file pointer and subRank");

  initialPos = 0;
  spacing    = 1;

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\n Finished computing filter parameters for chain " << m_name
                           << ": initialPos = " << initialPos
                           << ", spacing = "    << spacing
                           << "\n-----------------------------------------------------"
                           << "\n"
                           << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeHistKde( // Use the whole chain
  const uqChainStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                        passedOfs)
{
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n"
                           << "\n-----------------------------------------------------"
                           << "\n Computing histogram and/or KDE for chain " << m_name << " ..."
                           << "\n-----------------------------------------------------"
                           << "\n"
                           << std::endl;
  }

  int iRC = UQ_OK_RC;
  struct timeval timevalTmp;

  //****************************************************
  // Compute MIN and MAX: for histograms and KDE
  //****************************************************
  double tmpRunTime = 0.;
  iRC = gettimeofday(&timevalTmp, NULL);
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\nComputing min and max for histograms and KDE"
                           << std::endl;
  }

  V statsMinPositions(m_vectorSpace.zeroVector());
  V statsMaxPositions(m_vectorSpace.zeroVector());
  this->minMax(0, // Use the whole chain
               statsMinPositions,
               statsMaxPositions);

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\nComputed min values and max values for chain " << m_name
                           << std::endl;

    char line[512];
    sprintf(line,"%s",
            "Parameter");
    *m_env.subScreenFile() << line;

    sprintf(line,"%9s%s%9s%s",
            " ",
            "min",
            " ",
            "max");
    *m_env.subScreenFile() << line;

    for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
      sprintf(line,"\n%8.8s",
              m_vectorSpace.componentName(i).c_str() /*.*/);
      *m_env.subScreenFile() << line;

      sprintf(line,"%2s%11.4e%2s%11.4e",
              " ",
              statsMinPositions[i],
              " ",
              statsMaxPositions[i]);
      *m_env.subScreenFile() << line;
    }
    *m_env.subScreenFile() << std::endl;
  }

  V unifiedStatsMinPositions(statsMinPositions);
  V unifiedStatsMaxPositions(statsMaxPositions);
  if (m_env.numSubEnvironments() > 1) {
    // Compute unified min-max
    this->unifiedMinMax(0, // Use the whole chain
                        unifiedStatsMinPositions,
                        unifiedStatsMaxPositions);

    // Write unified min-max
    if (m_env.subScreenFile()) {
      if (m_vectorSpace.zeroVector().numberOfProcessorsRequiredForStorage() == 1) {
        if (m_env.inter0Rank() == 0) {
          *m_env.subScreenFile() << "\nComputed unified min values and max values for chain " << m_name
                                 << std::endl;

          char line[512];
          sprintf(line,"%s",
                  "Parameter");
          *m_env.subScreenFile() << line;

          sprintf(line,"%9s%s%9s%s",
                  " ",
                  "min",
                  " ",
                  "max");
          *m_env.subScreenFile() << line;

          for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
            sprintf(line,"\n%8.8s",
                    m_vectorSpace.componentName(i).c_str() /*.*/);
            *m_env.subScreenFile() << line;

            sprintf(line,"%2s%11.4e%2s%11.4e",
                    " ",
                    unifiedStatsMinPositions[i],
                    " ",
                    unifiedStatsMaxPositions[i]);
            *m_env.subScreenFile() << line;
          }
          *m_env.subScreenFile() << std::endl;
        }
      }
      else {
        UQ_FATAL_TEST_MACRO(true,
                            m_env.rank(),
                            "uqBaseVectorSequenceClass<V,M>::computeHistKde()",
                            "unified min-max writing, parallel vectors not supported yet");
      }
    } // if subScreenFile
  } // if numSubEnvs > 1

  m_env.fullComm().Barrier();
  tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "Chain min and max took " << tmpRunTime
                           << " seconds"
                           << std::endl;
  }

  //****************************************************
  // Compute histograms
  //****************************************************
  if ((statisticalOptions.histCompute()            ) &&
      (statisticalOptions.histNumInternalBins() > 0)) {
    tmpRunTime = 0.;
    iRC = gettimeofday(&timevalTmp, NULL);
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "\n-----------------------------------------------------"
                             << "\nComputing histograms"
                             << std::endl;
    }

    std::string subCoreName_HistCenters((std::string)(    "_HistCenters_sub")+m_env.subIdString());
    std::string uniCoreName_HistCenters((std::string)("_unifHistCenters_sub")+m_env.subIdString());
    if (m_env.numSubEnvironments() == 1) subCoreName_HistCenters = uniCoreName_HistCenters;

    std::string subCoreName_HistQuantts((std::string)(    "_HistQuantts_sub")+m_env.subIdString());
    std::string uniCoreName_HistQuantts((std::string)("_unifHistQuantts_sub")+m_env.subIdString());
    if (m_env.numSubEnvironments() == 1) subCoreName_HistQuantts = uniCoreName_HistQuantts;

    for (unsigned int i = 0; i < statsMaxPositions.size(); ++i) {
      statsMaxPositions[i] *= (1. + 1.e-15);
    }

    // Compute histograms
    std::vector<V*> histCentersForAllBins(statisticalOptions.histNumInternalBins()+2,NULL); // IMPORTANT: +2
    std::vector<V*> histQuanttsForAllBins(statisticalOptions.histNumInternalBins()+2,NULL); // IMPORTANT: +2
    this->histogram(0, // Use the whole chain
                    statsMinPositions,
                    statsMaxPositions,
                    histCentersForAllBins,
                    histQuanttsForAllBins);

    // Write histograms
    if (passedOfs) {
      std::ofstream& ofsvar = *passedOfs;
      ofsvar << m_name << subCoreName_HistCenters << " = zeros(" << this->vectorSize() /*.*/
             << ","                                              << histCentersForAllBins.size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int j = 0; j < histCentersForAllBins.size(); ++j) {
           ofsvar << m_name << subCoreName_HistCenters << "(" << i+1
                  << ","                                      << j+1
                  << ") = "                                   << (*(histCentersForAllBins[j]))[i]
                  << ";"
                  << std::endl;
        }
      }

      ofsvar << m_name << subCoreName_HistQuantts << " = zeros(" << this->vectorSize() /*.*/
             << ","                                              << histQuanttsForAllBins.size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int j = 0; j < histQuanttsForAllBins.size(); ++j) {
           ofsvar << m_name << subCoreName_HistQuantts << "(" << i+1
                  << ","                                      << j+1
                  << ") = "                                   << (*(histQuanttsForAllBins[j]))[i]
                  << ";"
                  << std::endl;
        }
      }
    }

    for (unsigned int i = 0; i < histQuanttsForAllBins.size(); ++i) {
      if (histQuanttsForAllBins[i] != NULL) delete histQuanttsForAllBins[i];
    }
    for (unsigned int i = 0; i < histCentersForAllBins.size(); ++i) {
      if (histCentersForAllBins[i] != NULL) delete histCentersForAllBins[i];
    }

    std::vector<V*> unifiedHistCentersForAllBins(statisticalOptions.histNumInternalBins()+2,NULL); // IMPORTANT: +2
    std::vector<V*> unifiedHistQuanttsForAllBins(statisticalOptions.histNumInternalBins()+2,NULL); // IMPORTANT: +2
    if (m_env.numSubEnvironments() > 1) {
      // Compute unified histogram
      this->unifiedHistogram(0, // Use the whole chain
                             unifiedStatsMinPositions,
                             unifiedStatsMaxPositions,
                             unifiedHistCentersForAllBins,
                             unifiedHistQuanttsForAllBins);

      // Write unified histogram
      if (passedOfs) {
        if (m_vectorSpace.zeroVector().numberOfProcessorsRequiredForStorage() == 1) {
          if (m_env.inter0Rank() == 0) {
            std::ofstream& ofsvar = *passedOfs;
            ofsvar << m_name << uniCoreName_HistCenters << " = zeros(" << this->vectorSize() /*.*/
                   << ","                                              << unifiedHistCentersForAllBins.size()
                   << ");"
                   << std::endl;
            for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
              for (unsigned int j = 0; j < unifiedHistCentersForAllBins.size(); ++j) {
                 ofsvar << m_name << uniCoreName_HistCenters << "(" << i+1
                        << ","                                      << j+1
                        << ") = "                                   << (*(unifiedHistCentersForAllBins[j]))[i]
                        << ";"
                        << std::endl;
              }
            }

            ofsvar << m_name << uniCoreName_HistQuantts << " = zeros(" << this->vectorSize() /*.*/
                   << ","                                              << unifiedHistQuanttsForAllBins.size()
                   << ");"
                   << std::endl;
            for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
              for (unsigned int j = 0; j < unifiedHistQuanttsForAllBins.size(); ++j) {
                 ofsvar << m_name << uniCoreName_HistQuantts << "(" << i+1
                        << ","                                      << j+1
                        << ") = "                                   << (*(unifiedHistQuanttsForAllBins[j]))[i]
                        << ";"
                        << std::endl;
              }
            }
          }
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_env.rank(),
                              "uqBaseVectorSequenceClass<V,M>::computeHistKde()",
                              "unified histogram writing, parallel vectors not supported yet");
        }
      } // if passedOfs

      for (unsigned int i = 0; i < unifiedHistQuanttsForAllBins.size(); ++i) {
        if (unifiedHistQuanttsForAllBins[i] != NULL) delete unifiedHistQuanttsForAllBins[i];
      }
      for (unsigned int i = 0; i < unifiedHistCentersForAllBins.size(); ++i) {
        if (unifiedHistCentersForAllBins[i] != NULL) delete unifiedHistCentersForAllBins[i];
      }
    } // if numSubEnvs > 1

    m_env.fullComm().Barrier();
    tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "Chain histograms took " << tmpRunTime
                             << " seconds"
                             << std::endl;
    }
  }

  //****************************************************
  // Compute estimations of probability densities
  //****************************************************
  if ((statisticalOptions.kdeCompute()             ) &&
      (statisticalOptions.kdeNumEvalPositions() > 0)) {
    tmpRunTime = 0.;
    iRC = gettimeofday(&timevalTmp, NULL);
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "\n-----------------------------------------------------"
                             << "\nComputing KDE"
                             << std::endl;
    }

    std::string subCoreName_GaussianKdePositions((std::string)(    "_GkdePosits_sub")+m_env.subIdString());
    std::string uniCoreName_GaussianKdePositions((std::string)("_unifGkdePosits_sub")+m_env.subIdString());
    if (m_env.numSubEnvironments() == 1) subCoreName_GaussianKdePositions = uniCoreName_GaussianKdePositions; // avoid temporarily (see '< -1' below)

    std::string subCoreName_GaussianKdeScaleVec ((std::string)(    "_GkdeScalev_sub")+m_env.subIdString());
    std::string uniCoreName_GaussianKdeScaleVec ((std::string)("_unifGkdeScalev_sub")+m_env.subIdString());
    if (m_env.numSubEnvironments() == 1) subCoreName_GaussianKdeScaleVec = uniCoreName_GaussianKdeScaleVec; // avoid temporarily (see '< -1' below)

    std::string subCoreName_GaussianKdeValues   ((std::string)(    "_GkdeValues_sub")+m_env.subIdString());
    std::string uniCoreName_GaussianKdeValues   ((std::string)("_unifGkdeValues_sub")+m_env.subIdString());
    if (m_env.numSubEnvironments() == 1) subCoreName_GaussianKdeValues = uniCoreName_GaussianKdeValues; // avoid temporarily (see '< -1' below)

    V iqrVec(m_vectorSpace.zeroVector());
    this->interQuantileRange(0, // Use the whole chain
                             iqrVec);

    V gaussianKdeScaleVec(m_vectorSpace.zeroVector());
    this->scalesForKDE(0, // Use the whole chain
                       iqrVec,
                       gaussianKdeScaleVec);

    std::vector<V*> kdeEvalPositions(statisticalOptions.kdeNumEvalPositions(),NULL);
    uqMiscComputePositionsBetweenMinMax(statsMinPositions,
                                        statsMaxPositions,
                                        kdeEvalPositions);

    std::vector<V*> gaussianKdeDensities(statisticalOptions.kdeNumEvalPositions(),NULL);
    this->gaussianKDE(0, // Use the whole chain
                      gaussianKdeScaleVec,
                      kdeEvalPositions,
                      gaussianKdeDensities);

    // Write iqr
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "\nComputed inter quantile ranges for chain " << m_name
                               << std::endl;

      char line[512];
      sprintf(line,"%s",
              "Parameter");
      *m_env.subScreenFile() << line;

      sprintf(line,"%9s%s",
              " ",
              "iqr");
      *m_env.subScreenFile() << line;

      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        sprintf(line,"\n%8.8s",
                m_vectorSpace.componentName(i).c_str() /*.*/);
        *m_env.subScreenFile() << line;

        sprintf(line,"%2s%11.4e",
                " ",
                iqrVec[i]);
        *m_env.subScreenFile() << line;
      }
      *m_env.subScreenFile() << std::endl;
    }

    // Write KDE
    if (passedOfs) {
      std::ofstream& ofsvar = *passedOfs;
      ofsvar << m_name << subCoreName_GaussianKdePositions << " = zeros(" << this->vectorSize() /*.*/
             << ","                                                       << kdeEvalPositions.size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int j = 0; j < kdeEvalPositions.size(); ++j) {
          ofsvar << m_name << subCoreName_GaussianKdePositions << "(" << i+1
                 << ","                                               << j+1
                 << ") = "                                            << (*(kdeEvalPositions[j]))[i]
                 << ";"
                 << std::endl;
        }
      }

      ofsvar << m_name << subCoreName_GaussianKdeScaleVec << " = zeros(" << this->vectorSize() /*.*/
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        ofsvar << m_name << subCoreName_GaussianKdeScaleVec << "(" << i+1
               << ") = "                                           << gaussianKdeScaleVec[i]
               << ";"
               << std::endl;
      }

      ofsvar << m_name << subCoreName_GaussianKdeValues << " = zeros(" << this->vectorSize() /*.*/
             << ","                                                    << gaussianKdeDensities.size()
             << ");"
             << std::endl;
      for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
        for (unsigned int j = 0; j < gaussianKdeDensities.size(); ++j) {
          ofsvar << m_name << subCoreName_GaussianKdeValues << "(" << i+1
                 << ","                                            << j+1
                 << ") = "                                         << (*(gaussianKdeDensities[j]))[i]
                 << ";"
                 << std::endl;
        }
      }
    }

    for (unsigned int i = 0; i < gaussianKdeDensities.size(); ++i) {
      if (gaussianKdeDensities[i] != NULL) delete gaussianKdeDensities[i];
    }
    for (unsigned int i = 0; i < kdeEvalPositions.size(); ++i) {
      if (kdeEvalPositions[i] != NULL) delete kdeEvalPositions[i];
    }

    if ((int) m_env.numSubEnvironments() > 1) { // avoid code temporarily
      // Compute unified KDE
      V unifiedIqrVec(m_vectorSpace.zeroVector());
      this->unifiedInterQuantileRange(0, // Use the whole chain
                                      unifiedIqrVec);
      m_env.fullComm().Barrier();

      V unifiedGaussianKdeScaleVec(m_vectorSpace.zeroVector());
      this->unifiedScalesForKDE(0, // Use the whole chain
                                unifiedIqrVec,
                                unifiedGaussianKdeScaleVec);
      m_env.fullComm().Barrier();

      std::vector<V*> unifiedKdeEvalPositions(statisticalOptions.kdeNumEvalPositions(),NULL);
      uqMiscComputePositionsBetweenMinMax(unifiedStatsMinPositions,
                                          unifiedStatsMaxPositions,
                                          unifiedKdeEvalPositions);

      std::vector<V*> unifiedGaussianKdeDensities(statisticalOptions.kdeNumEvalPositions(),NULL);
      this->unifiedGaussianKDE(0, // Use the whole chain
                               unifiedGaussianKdeScaleVec,
                               unifiedKdeEvalPositions,
                               unifiedGaussianKdeDensities);
      m_env.fullComm().Barrier();

      // Write unified iqr
      if (m_env.subScreenFile()) {
        if (m_vectorSpace.zeroVector().numberOfProcessorsRequiredForStorage() == 1) {
          if (m_env.inter0Rank() == 0) {
            *m_env.subScreenFile() << "\nComputed unified inter quantile ranges for chain " << m_name
                                   << std::endl;

            char line[512];
            sprintf(line,"%s",
                    "Parameter");
            *m_env.subScreenFile() << line;

            sprintf(line,"%9s%s",
                    " ",
                    "iqr");
            *m_env.subScreenFile() << line;

            for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
              sprintf(line,"\n%8.8s",
                      m_vectorSpace.componentName(i).c_str() /*.*/);
              *m_env.subScreenFile() << line;

              sprintf(line,"%2s%11.4e",
                      " ",
                      unifiedIqrVec[i]);
              *m_env.subScreenFile() << line;
            }
            *m_env.subScreenFile() << std::endl;
          }
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_env.rank(),
                              "uqBaseVectorSequenceClass<V,M>::computeHistKde()",
                              "unified iqr writing, parallel vectors not supported yet");
        }
      }

      // Write unified KDE
      if (passedOfs) {
        if (m_vectorSpace.zeroVector().numberOfProcessorsRequiredForStorage() == 1) {
          if (m_env.inter0Rank() == 0) {
            std::ofstream& ofsvar = *passedOfs;
            ofsvar << m_name << uniCoreName_GaussianKdePositions << " = zeros(" << this->vectorSize() /*.*/
                   << ","                                                       << unifiedKdeEvalPositions.size()
                   << ");"
                   << std::endl;
            for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
              for (unsigned int j = 0; j < unifiedKdeEvalPositions.size(); ++j) {
                ofsvar << m_name << uniCoreName_GaussianKdePositions << "(" << i+1
                       << ","                                               << j+1
                       << ") = "                                            << (*(unifiedKdeEvalPositions[j]))[i]
                       << ";"
                       << std::endl;
              }
            }

            ofsvar << m_name << uniCoreName_GaussianKdeScaleVec << " = zeros(" << this->vectorSize() /*.*/
                   << ");"
                   << std::endl;
            for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
              ofsvar << m_name << uniCoreName_GaussianKdeScaleVec << "(" << i+1
                     << ") = "                                           << unifiedGaussianKdeScaleVec[i]
                     << ";"
                     << std::endl;
            }

            ofsvar << m_name << uniCoreName_GaussianKdeValues << " = zeros(" << this->vectorSize() /*.*/
                   << ","                                                    << unifiedGaussianKdeDensities.size()
                   << ");"
                   << std::endl;
            for (unsigned int i = 0; i < this->vectorSize() /*.*/; ++i) {
              for (unsigned int j = 0; j < unifiedGaussianKdeDensities.size(); ++j) {
                ofsvar << m_name << uniCoreName_GaussianKdeValues << "(" << i+1
                       << ","                                            << j+1
                       << ") = "                                         << (*(unifiedGaussianKdeDensities[j]))[i]
                       << ";"
                       << std::endl;
              }
            }
          }
        }
        else {
          UQ_FATAL_TEST_MACRO(true,
                              m_env.rank(),
                              "uqBaseVectorSequenceClass<V,M>::computeHistKde()",
                              "unified KDE writing, parallel vectors not supported yet");
        }
      }

      for (unsigned int i = 0; i < unifiedGaussianKdeDensities.size(); ++i) {
        if (unifiedGaussianKdeDensities[i] != NULL) delete unifiedGaussianKdeDensities[i];
      }
      for (unsigned int i = 0; i < unifiedKdeEvalPositions.size(); ++i) {
        if (unifiedKdeEvalPositions[i] != NULL) delete unifiedKdeEvalPositions[i];
      }
    }

    m_env.fullComm().Barrier();
    tmpRunTime += uqMiscGetEllapsedSeconds(&timevalTmp);
    if (m_env.subScreenFile()) {
      *m_env.subScreenFile() << "Chain KDE took " << tmpRunTime
                             << " seconds"
                             << std::endl;
    }
  }

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\n Finished computing histogram and/or KDE for chain " << m_name
                           << "\n-----------------------------------------------------"
                           << "\n"
                           << std::endl;
  }

  return;
}

template<class V, class M>
void
uqBaseVectorSequenceClass<V,M>::computeCovCorrMatrices( // Use the whole chain
  const uqChainStatisticalOptionsClass& statisticalOptions,
  std::ofstream*                        passedOfs)
{
  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n"
                           << "\n-----------------------------------------------------"
                           << "\n Computing covariance and correlation matrices for chain " << m_name << " ..."
                           << "\n-----------------------------------------------------"
                           << "\n"
                           << std::endl;
  }

  //int iRC = UQ_OK_RC;
  //struct timeval timevalTmp;
  M* covarianceMatrix = new M(m_env,
                              m_vectorSpace.map(),  // number of rows
                              m_vectorSpace.dim()); // number of cols
  M* correlationMatrix = new M(m_env,
                               m_vectorSpace.map(),  // number of rows
                               m_vectorSpace.dim()); // number of cols

  uqComputeCovCorrMatricesBetweenVectorSequences(*this,
                                                 *this,
                                                 this->sequenceSize(),
                                                 *covarianceMatrix,
                                                 *correlationMatrix);

  if (m_env.subScreenFile()) {
    if (m_vectorSpace.zeroVector().numberOfProcessorsRequiredForStorage() == 1) {
      // Only unified covariance matrix is written. And only one processor writes it.
      if (m_env.inter0Rank() == 0) {
        *m_env.subScreenFile() << "\nuqBaseVectorSequenceClass<V,M>::computeCovCorrMatrices"
                               << ", chain " << m_name
                               << ": contents of covariance matrix are\n" << *covarianceMatrix
                               << std::endl;

        *m_env.subScreenFile() << "\nuqBaseVectorSequenceClass<V,M>::computeCovCorrMatrices"
                               << ", chain " << m_name
                               << ": contents of correlation matrix are\n" << *correlationMatrix
                               << std::endl;
      }
    }
    else {
      UQ_FATAL_TEST_MACRO(true,
                          m_env.rank(),
                          "uqBaseVectorSequenceClass<V,M>::computeCovCorrMatrices()",
                          "parallel vectors not supported yet");
    }
  }

  if (m_env.subScreenFile()) {
    *m_env.subScreenFile() << "\n-----------------------------------------------------"
                           << "\n Finished computing covariance and correlation matrices for chain " << m_name
                           << "\n-----------------------------------------------------"
                           << "\n"
                           << std::endl;
  }

  return;
}

template <class P_V, class P_M, class Q_V, class Q_M>
void
uqComputeCovCorrMatricesBetweenVectorSequences(
  const uqBaseVectorSequenceClass<P_V,P_M>& localPSeq,
  const uqBaseVectorSequenceClass<Q_V,Q_M>& localQSeq,
        unsigned int                        localNumSamples,
        P_M&                                pqCovMatrix,
        P_M&                                pqCorrMatrix)
{
  // Check input data consistency
  const uqBaseEnvironmentClass& env = localPSeq.vectorSpace().zeroVector().env();

  bool useOnlyInter0Comm = (localPSeq.vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage() == 1) &&
                           (localQSeq.vectorSpace().zeroVector().numberOfProcessorsRequiredForStorage() == 1);

  UQ_FATAL_TEST_MACRO((useOnlyInter0Comm == false),
                      env.rank(),
                      "uqComputeCovCorrMatricesBetweenVectorSequences()",
                      "parallel vectors not supported yet");

  unsigned int numRows = localPSeq.vectorSpace().dim();
  unsigned int numCols = localQSeq.vectorSpace().dim();

  UQ_FATAL_TEST_MACRO((numRows != pqCovMatrix.numRows()) || (numCols != pqCovMatrix.numCols()),
                      env.rank(),
                      "uqComputeCovCorrMatricesBetweenVectorSequences()",
                      "inconsistent dimensions for covariance matrix");

  UQ_FATAL_TEST_MACRO((numRows != pqCorrMatrix.numRows()) || (numCols != pqCorrMatrix.numCols()),
                      env.rank(),
                      "uqComputeCorrelationBetweenVectorSequences()",
                      "inconsistent dimensions for correlation matrix");

  UQ_FATAL_TEST_MACRO((localNumSamples > localPSeq.sequenceSize()) || (localNumSamples > localQSeq.sequenceSize()),
                      env.rank(),
                      "uqComputeCovCorrMatricesBetweenVectorSequences()",
                      "localNumSamples is too large");

  // For both P and Q vector sequences: fill them
  P_V tmpP(localPSeq.vectorSpace().zeroVector());
  Q_V tmpQ(localQSeq.vectorSpace().zeroVector());

  // For both P and Q vector sequences: compute the unified mean
  P_V unifiedMeanP(localPSeq.vectorSpace().zeroVector());
  localPSeq.unifiedMean(0,localNumSamples,unifiedMeanP);

  Q_V unifiedMeanQ(localQSeq.vectorSpace().zeroVector());
  localQSeq.unifiedMean(0,localNumSamples,unifiedMeanQ);

  // Compute "local" covariance matrix
  for (unsigned i = 0; i < numRows; ++i) {
    for (unsigned j = 0; j < numCols; ++j) {
      pqCovMatrix(i,j) = 0.;
    }
  }
  for (unsigned k = 0; k < localNumSamples; ++k) {
    // For both P and Q vector sequences: get the difference (wrt the unified mean) in them
    localPSeq.getPositionValues(k,tmpP);
    tmpP -= unifiedMeanP;

    localQSeq.getPositionValues(k,tmpQ);
    tmpQ -= unifiedMeanQ;

    for (unsigned i = 0; i < numRows; ++i) {
      for (unsigned j = 0; j < numCols; ++j) {
        pqCovMatrix(i,j) += tmpP[i]*tmpQ[j];
      }
    }
  }

  // For both P and Q vector sequences: compute the unified variance
  P_V unifiedSampleVarianceP(localPSeq.vectorSpace().zeroVector());
  localPSeq.unifiedSampleVariance(0,
                                  localNumSamples,
                                  unifiedMeanP,
                                  unifiedSampleVarianceP);

  Q_V unifiedSampleVarianceQ(localQSeq.vectorSpace().zeroVector());
  localQSeq.unifiedSampleVariance(0,
                                  localNumSamples,
                                  unifiedMeanQ,
                                  unifiedSampleVarianceQ);

  // Compute unified covariance matrix
  if (useOnlyInter0Comm) {
    if (env.inter0Rank() >= 0) {
      unsigned int unifiedNumSamples = 0;
      int mpiRC = MPI_Allreduce((void *) &localNumSamples, (void *) &unifiedNumSamples, (int) 1, MPI_UNSIGNED, MPI_SUM, env.inter0Comm().Comm());
      UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                          env.rank(),
                          "uqComputeCovCorrMatricesBetweenVectorSequences()",
                          "failed MPI_Allreduce() for localNumSamples");

      for (unsigned i = 0; i < numRows; ++i) {
        for (unsigned j = 0; j < numCols; ++j) {
          double aux = 0.;
          int mpiRC = MPI_Allreduce((void *) &pqCovMatrix(i,j), (void *) &aux, (int) 1, MPI_DOUBLE, MPI_SUM, env.inter0Comm().Comm());
          UQ_FATAL_TEST_MACRO(mpiRC != MPI_SUCCESS,
                              env.rank(),
                              "uqComputeCovCorrMatricesBetweenVectorSequences()",
                              "failed MPI_Allreduce() for a matrix position");
          pqCovMatrix(i,j) = aux/((double) unifiedNumSamples);
        }
      }

      for (unsigned i = 0; i < numRows; ++i) {
        for (unsigned j = 0; j < numCols; ++j) {
          pqCorrMatrix(i,j) = pqCovMatrix(i,j)/sqrt(unifiedSampleVarianceP[i])/sqrt(unifiedSampleVarianceQ[j]);
        }
      }
    }
    else {
      // Node not in the 'inter0' communicator: do nothing extra
    }
  }
  else {
    UQ_FATAL_TEST_MACRO((useOnlyInter0Comm == false),
                        env.rank(),
                        "uqComputeCovCorrMatricesBetweenVectorSequences()",
                        "parallel vectors not supported yet (2)");
  }

  return;
}
#endif // __UQ_VECTOR_SEQUENCE_H__

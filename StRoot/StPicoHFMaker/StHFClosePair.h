#ifndef StHFClosePair_H
#define StHFClosePair_H

#include "TObject.h"
#include "TClonesArray.h"

#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StarClassLibrary/StThreeVectorF.hh"

class StPicoTrack;

class StHFClosePair : public TObject
{
public:
  StHFClosePair();
  StHFClosePair(StHFClosePair const &);
  StHFClosePair(StPicoTrack const * particle1, StPicoTrack const * particle2, 
		unsigned short p1Idx, unsigned short p2Idx,
		StThreeVectorF const & vtx, float bField, bool useStraightLine = true);

  StHFClosePair & operator= (StHFClosePair const &);
  ~StHFClosePair() {;}

  // getters
  float particle1Dca() const;
  float particle2Dca() const;
  float dcaDaughters() const;
  StThreeVectorF p1AtDcaToP2() const;
  StThreeVectorF p2AtDcaToP1() const;
  unsigned short particle1Idx() const;
  unsigned short particle2Idx() const;  
  float p1massHypothesis() const;
  float p2massHypothesis() const;

private:
  float mParticle1Dca;
  float mParticle2Dca;
  float mDcaDaughters;

  StThreeVectorF mP1AtDcaToP2;
  StThreeVectorF mP2AtDcaToP1;

  float mMassHypothesis1;
  float mMassHypothesis2;

  unsigned short  mParticle1Idx; // index of track in StPicoDstEvent 
  unsigned short  mParticle2Idx; // index of track in StPicoDstEvent for particle, idx in tertiary vertex array for pair 
};

// getters
inline float StHFClosePair::particle1Dca() const { return mParticle1Dca; }
inline float StHFClosePair::particle2Dca() const { return mParticle2Dca; }
inline float StHFClosePair::StdcaDaughters() const { return mDcaDaughters; }
inline StThreeVectorF StHFClosePair::p1AtDcaToP2() const { return mP1AtDcaToP2; }
inline StThreeVectorF StHFClosePair::p2AtDcaToP1() const { return mP2AtDcaToP1; }
inline unsigned short StHFClosePair::particle1Idx() const { return mParticle1Dca; }
inline unsigned short StHFClosePair::particle2Idx() const { return mParticle2Dca; }
inline float StHFClosePair::p1massHypothesis() const { return mMassHypothesis1; }
inline float StHFClosePair::p2massHypothesis() const { return mMassHypothesis2; }

#endif

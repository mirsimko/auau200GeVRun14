#ifndef StHFClosePair_HH
#define StHFClosePair_HH

/* **************************************************
 *  Class based on StPair created for three body decays 
 *  
 *  The purpose of this class, compared to StHFPair, is to not compute things that 
 *  aren't necessary for the three body decays like mass, etc. Objects like straight
 *  line approximation of the tracks are kept for later use for the StHFTriplet
 * **************************************************
 *
 *  Author: 
 *        **  Miroslav Simko
 *
 *  ** Code Maintainer 
 *
 * **************************************************
 */
#include "TObject.h"
#include "TClonesArray.h"

#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"

class StPicoTrack;

class StHFClosePair : public TObject
{
public:
  StHFClosePair();
  StHFClosePair(StHFClosePair const &);
  StHFClosePair(StPicoTrack const * particle1, StPicoTrack const * particle2, 
		float p1mass, float p2mass,
		unsigned short p1Idx, unsigned short p2Idx,
		StThreeVectorF const & vtx, float bField, bool useStraightLine = true);

  StHFClosePair & operator= (StHFClosePair const &rhs)
  {
    StHFClosePair tmp(rhs);
    std::swap(*this, tmp);
    return *this;
  }
  ~StHFClosePair() ;

  // getters
  float particle1Dca() const;
  float particle2Dca() const;
  float dcaDaughters() const;
  StThreeVectorF p1AtDcaToP2() const;
  StThreeVectorF p2AtDcaToP1() const;
  StPhysicalHelixD * p1StraightLine();
  StPhysicalHelixD * p2StraightLine();
  unsigned short particle1Idx() const;
  unsigned short particle2Idx() const;  
  float p1massHypothesis() const;
  float p2massHypothesis() const;
  StPhysicalHelixD * p1Helix();
  StPhysicalHelixD * p2Helix();

private:
  float mParticle1Dca;
  float mParticle2Dca;
  float mDcaDaughters;

  unsigned short  mParticle1Idx; // index of track in StPicoDstEvent 
  unsigned short  mParticle2Idx; // index of track in StPicoDstEvent for particle, idx in tertiary vertex array for pair 
  float mMassHypothesis1;
  float mMassHypothesis2;

  StThreeVectorF mP1AtDcaToP2;
  StThreeVectorF mP2AtDcaToP1;
  StPhysicalHelixD * mP1StraightLine;
  StPhysicalHelixD * mP2StraightLine;
  StPhysicalHelixD * mP1Helix ;
  StPhysicalHelixD * mP2Helix ;

  ClassDef(StHFClosePair, 2)
};

// getters
inline float StHFClosePair::particle1Dca() const { return mParticle1Dca; }
inline float StHFClosePair::particle2Dca() const { return mParticle2Dca; }
inline float StHFClosePair::dcaDaughters() const { return mDcaDaughters; }
inline unsigned short StHFClosePair::particle1Idx() const { return mParticle1Dca; }
inline unsigned short StHFClosePair::particle2Idx() const { return mParticle2Dca; }
inline float StHFClosePair::p1massHypothesis() const { return mMassHypothesis1; }
inline float StHFClosePair::p2massHypothesis() const { return mMassHypothesis2; }
inline StThreeVectorF StHFClosePair::p1AtDcaToP2() const { return mP1AtDcaToP2; }
inline StThreeVectorF StHFClosePair::p2AtDcaToP1() const { return mP2AtDcaToP1; }
inline StPhysicalHelixD * StHFClosePair::p1StraightLine() { return mP1StraightLine; }
inline StPhysicalHelixD * StHFClosePair::p2StraightLine() { return mP2StraightLine; }
inline StPhysicalHelixD * StHFClosePair::p1Helix() { return mP1Helix; }
inline StPhysicalHelixD * StHFClosePair::p2Helix() { return mP2Helix; }

#endif

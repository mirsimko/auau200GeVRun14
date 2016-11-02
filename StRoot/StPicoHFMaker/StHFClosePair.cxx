#include <limits>
#include <cmath>
#include <utility>

#include "StHFPair.h"

#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StPicoDstMaker/StPicoTrack.h"

ClassImp(StHFPair)

// _________________________________________________________
StHFClosePair::StHFClosePair() : mDecayLength(std::numeric_limits<float>::quiet_NaN()), 
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()), 
  mDcaDaughters(std::numeric_limits<float>::max()), 
  mParticle1Idx(std::numeric_limits<unsigned short>::max()), mParticle2Idx(std::numeric_limits<unsigned short>::max()), 
  mMassHypothesis1(std::numeric_limits<float>::quiet_NaN()), mMassHypothesis2(std::numeric_limits<float>::quiet_NaN())
{}

// _________________________________________________________
StHFClosePair::StHFClosePair(StHFClosePair const &rhs) : mDecayLength(rhs.decayLength),
  mParticle1Dca(rhs.mParticle1Dca), mParticle2Dca(rhs.mParticle2Dca), mDcaToPrimaryVertex(rhs.mDcaToPrimaryVertex),
  mDcaDaughters(rhs.mDcaDaughters), mParticle1Idx(rhs.mParticle1Idx), mParticle2Idx(rhs.mParticle2Idx),
  mMassHypothesis1(rhs.mMassHypothesis1), mMassHypothesis2(rhs.mMassHypothesis2),
  mP1AtDcaToP2(rhs.mP1AtDcaToP2), mP2AtDcaToP1(rhs.mP2AtDcaToP1)
{}

// _________________________________________________________
StHFClosePair & operator= (StHFClosePair const &rhs)
{
  StHFClosePair tmp(rhs);
  std::swap(*this, tmp);
  return *this;
}

// _________________________________________________________
StHFClosePair::StHFClosePair(StPicoTrack const * particle1, StPicoTrack const * particle2, 
			     unsigned short p1Idx, unsigned short p2Idx,
			     StThreeVectorF const & vtx, float bField, bool useStraightLine = true) :
  mDecayLength(std::numeric_limits<float>::quiet_NaN()), 
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()), 
  mDcaToPrimaryVertex(std::numeric_limits<float>::quiet_NaN()), mDcaDaughters(std::numeric_limits<float>::max()), 
  mParticle1Idx(p1Idx), mParticle2Idx(p2Idx)
{
  // see if the particles are the same
  if ((!particle1 || !particle2) || (mParticle1Idx == mParticle2Idx) || (particle1->id() == particle2->id())) {
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    return;
  }
  // -- move origins of helices to the primary vertex origin
  p1Helix.moveOrigin(p1Helix.pathLength(vtx));
  p2Helix.moveOrigin(p2Helix.pathLength(vtx));

  // -- use straight lines approximation to get point of DCA of particle1-particle2 pair
  StThreeVectorF const p1Mom = p1Helix.momentum(bField * kilogauss);
  StThreeVectorF const p2Mom = p2Helix.momentum(bField * kilogauss);
  StPhysicalHelixD const p1StraightLine(p1Mom, p1Helix.origin(), 0, particle1->charge());
  StPhysicalHelixD const p2StraightLine(p2Mom, p2Helix.origin(), 0, particle2->charge());

  pair<double, double> const ss = (useStraightLine) ? p1StraightLine.pathLengths(p2StraightLine) : p1Helix.pathLengths(p2Helix);
  mP1AtDcaToP2 = (useStraightLine) ? p1StraightLine.at(ss.first) : p1Helix.at(ss.first);
  mP2AtDcaToP1 = (useStraightLine) ? p2StraightLine.at(ss.second) : p2Helix.at(ss.second);

  // -- calculate DCA of particle1 to particle2 at their DCA
  mDcaDaughters = (mP1AtDcaToP2 - mP2AtDcaToP1).mag();

  // -- single part DCAs
  mParticle1Dca = (p1Helix.origin() - vtx).mag();
  mParticle2Dca = (p2Helix.origin() - vtx).mag();
}

// _________________________________________________________

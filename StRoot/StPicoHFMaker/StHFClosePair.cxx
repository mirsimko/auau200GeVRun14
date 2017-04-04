#include <limits>
#include <cmath>
#include <utility>
#include <iostream>

#include "StHFClosePair.h"

#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StPicoDstMaker/StPicoTrack.h"

ClassImp(StHFClosePair)

// _________________________________________________________
StHFClosePair::StHFClosePair() :
  mP1Helix(NULL), 
  mP2Helix(NULL),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), 
  mParticle2Dca(std::numeric_limits<float>::quiet_NaN()), 
  mDcaDaughters(std::numeric_limits<float>::max()), 
  mParticle1Idx(std::numeric_limits<unsigned short>::max()), 
  mParticle2Idx(std::numeric_limits<unsigned short>::max()), 
  mMassHypothesis1(std::numeric_limits<float>::quiet_NaN()), 
  mMassHypothesis2(std::numeric_limits<float>::quiet_NaN()),
  mP1StraightLine(NULL), 
  mP2StraightLine(NULL) 
{}

// _________________________________________________________
StHFClosePair::StHFClosePair(StHFClosePair const &rhs) :
  mParticle1Dca(rhs.mParticle1Dca), 
  mParticle2Dca(rhs.mParticle2Dca),
  mDcaDaughters(rhs.mDcaDaughters), 
  mParticle1Idx(rhs.mParticle1Idx), 
  mParticle2Idx(rhs.mParticle2Idx),
  mMassHypothesis1(rhs.mMassHypothesis1), 
  mMassHypothesis2(rhs.mMassHypothesis2),
  mP1AtDcaToP2(rhs.mP1AtDcaToP2), 
  mP2AtDcaToP1(rhs.mP2AtDcaToP1)
{
  mP1StraightLine = new StPhysicalHelixD(*(rhs.mP1StraightLine)); // do a deep copy of the pointees
  mP2StraightLine = new StPhysicalHelixD(*(rhs.mP2StraightLine));
  mP1Helix = new StPhysicalHelixD(*(rhs.mP1Helix));
  mP2Helix = new StPhysicalHelixD(*(rhs.mP2Helix));
}

// _________________________________________________________
StHFClosePair::~StHFClosePair()
{
  delete mP1StraightLine;
  delete mP2StraightLine;
  delete mP1Helix;
  delete mP2Helix;
}

// _________________________________________________________
StHFClosePair::StHFClosePair(StPicoTrack const * particle1, StPicoTrack const * particle2, 
			     float p1mass, float p2mass,
			     unsigned short p1Idx, unsigned short p2Idx,
			     StThreeVectorF const & vtx, float bField, bool useStraightLine) :
  mP1Helix(NULL), 
  mP2Helix(NULL),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), 
  mParticle2Dca(std::numeric_limits<float>::quiet_NaN()), 
  mDcaDaughters(std::numeric_limits<float>::max()), 
  mMassHypothesis1(p1mass), 
  mMassHypothesis2(p2mass),
  mParticle1Idx(p1Idx), 
  mParticle2Idx(p2Idx),
  mP1StraightLine(NULL), 
  mP2StraightLine(NULL)
{
  // see if the particles are the same
  if(!particle1 || !particle2 || mParticle1Idx == mParticle2Idx ||  particle1->id() == particle2->id()) {
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    return;
  }

  mP1Helix = new StPhysicalHelixD( particle1->dcaGeometry().helix());
  mP2Helix = new StPhysicalHelixD( particle2->dcaGeometry().helix());
  if (!mP1Helix || !mP2Helix)
  {
    cerr << "StHFClosePair::StHFClosePair(...): Helices not initiated" << endl;
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    return;
  }
  calculateTopology(mP1Helix, mP2Helix, p1mass, p2mass, p1Idx, p2Idx, vtx, useStraightLine);
}

// _________________________________________________________
void StHFClosePair::calculateTopology(StPhysicalHelixD *p1Helix, StPhysicalHelixD *p2Helix, 
				      float p1mass, float p2mass,
				      unsigned short p1Idx, unsigned short p2Idx,
				      StThreeVectorF const & vtx, float bField, bool useStraightLine = true):
  mP1Helix(p1Helix), 
  mP2Helix(p2Helix),
  mMassHypothesis1(p1mass), 
  mMassHypothesis2(p2mass),
  mParticle1Idx(p1Idx), 
  mParticle2Idx(p2Idx)
{
  if (!mP1Helix || !mP2Helix)
  {
    cerr << "StHFClosePair::calculateTopology(...): Helices not found" << endl;
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    return;
  }
  // -- move origins of helices to the primary vertex origin
  mP1Helix->moveOrigin(mP1Helix->pathLength(vtx));
  mP2Helix->moveOrigin(mP2Helix->pathLength(vtx));

  // -- use straight lines approximation to get point of DCA of particle1-particle2 pair
  StThreeVectorF const p1Mom = mP1Helix->momentum(bField * kilogauss);
  StThreeVectorF const p2Mom = mP2Helix->momentum(bField * kilogauss);
  mP1StraightLine = new StPhysicalHelixD(p1Mom, mP1Helix->origin(), 0, particle1->charge());
  mP2StraightLine = new StPhysicalHelixD(p2Mom, mP2Helix->origin(), 0, particle2->charge());

  if(!mP1StraightLine || !mP2StraightLine)
  {
    cerr << "StHFClosePair::StHFClosePair(...): StraightLines not initiated" << endl;
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    return;
  }
  pair<double, double> const ss = (useStraightLine) ? mP1StraightLine->pathLengths(*mP2StraightLine) : mP1Helix->pathLengths(*mP2Helix);
  mP1AtDcaToP2 = (useStraightLine) ? mP1StraightLine->at(ss.first) : mP1Helix->at(ss.first);
  mP2AtDcaToP1 = (useStraightLine) ? mP2StraightLine->at(ss.second) : mP2Helix->at(ss.second);

  // -- calculate DCA of particle1 to particle2 at their DCA
  mDcaDaughters = (mP1AtDcaToP2 - mP2AtDcaToP1).mag();

  // -- single part DCAs
  mParticle1Dca = (mP1Helix->origin() - vtx).mag();
  mParticle2Dca = (mP2Helix->origin() - vtx).mag();
}

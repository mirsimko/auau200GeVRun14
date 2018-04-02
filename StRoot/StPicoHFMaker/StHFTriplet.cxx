#include <limits>
#include <cmath>

#include "StHFTriplet.h"

#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoHFMaker/StHFClosePair.h"

ClassImp(StHFTriplet)

// _________________________________________________________
StHFTriplet::StHFTriplet(): mLorentzVector(StLorentzVectorF()), mDecayVertex(StThreeVectorF()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()), 
  mParticle3Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(std::numeric_limits<unsigned short>::max()), mParticle2Idx(std::numeric_limits<unsigned short>::max()), 
  mParticle3Idx(std::numeric_limits<unsigned short>::max()),
  mDcaDaughters12(std::numeric_limits<float>::max()),  mDcaDaughters23(std::numeric_limits<float>::max()),
  mDcaDaughters31(std::numeric_limits<float>::max()), mDV0Max(std::numeric_limits<float>::quiet_NaN()){
}

// _________________________________________________________
StHFTriplet::StHFTriplet(StHFTriplet const * t) : 
  mLorentzVector(t->mLorentzVector), mDecayVertex(t->mDecayVertex),
  mPointingAngle(t->mPointingAngle), mDecayLength(t->mDecayLength), 
  mParticle1Dca(t->mParticle1Dca), mParticle2Dca(t->mParticle2Dca), mParticle3Dca(t->mParticle3Dca),
  mParticle1Idx(t->mParticle1Idx), mParticle2Idx(t->mParticle2Idx), mParticle3Idx(t->mParticle3Idx),
  mDcaDaughters12(t->mDcaDaughters12),  mDcaDaughters23(t->mDcaDaughters23), mDcaDaughters31(t->mDcaDaughters31),
  mDV0Max(t->mDV0Max){
}

// _________________________________________________________
StHFTriplet::StHFTriplet(StPicoTrack const * const particle1, StPicoTrack const * const particle2, StPicoTrack const * const particle3,
			 float p1MassHypo, float p2MassHypo, float p3MassHypo,
			 unsigned short const p1Idx, unsigned short const p2Idx, unsigned short const p3Idx,
			 StThreeVectorF const & vtx, float const bField)  : 
  mLorentzVector(StLorentzVectorF()), mDecayVertex(StThreeVectorF()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(std::numeric_limits<float>::quiet_NaN()), mParticle2Dca(std::numeric_limits<float>::quiet_NaN()), 
  mParticle3Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(p1Idx), mParticle2Idx(p2Idx),  mParticle3Idx(p3Idx),
  mDcaDaughters12(std::numeric_limits<float>::max()), mDcaDaughters23(std::numeric_limits<float>::max()),  
  mDcaDaughters31(std::numeric_limits<float>::max()),
  mDV0Max(std::numeric_limits<float>::quiet_NaN())
{
  // implemented via StHFClosePair
  if ((!particle1 || !particle2 || !particle3) || 
      (particle1->id() == particle2->id() || particle1->id() == particle3->id() || particle2->id() == particle3->id())) {
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    mParticle3Idx = std::numeric_limits<unsigned short>::max();
    return;
  }

  StHFClosePair closePair(particle1, particle2, p1MassHypo, p2MassHypo, p1Idx, p2Idx, vtx, bField);

  StPhysicalHelixD p3Helix = particle3->helix(bField);

  calculateTopology(&closePair, p3Helix, p3MassHypo, particle3->charge(), p3Idx, vtx, bField);
}

// _________________________________________________________
StHFTriplet::StHFTriplet(StHFClosePair * closePair, StPicoTrack const * particle3, 
			float p3MassHypo,
			unsigned short p3Idx,
			StThreeVectorF const & vtx, float bField) :
  mLorentzVector(StLorentzVectorF()), 
  mDecayVertex(StThreeVectorF()),
  mPointingAngle(std::numeric_limits<float>::quiet_NaN()), 
  mDecayLength(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Dca(closePair->particle1Dca()), 
  mParticle2Dca(closePair->particle2Dca()), 
  mParticle3Dca(std::numeric_limits<float>::quiet_NaN()),
  mParticle1Idx(closePair->particle1Idx()), 
  mParticle2Idx(closePair->particle2Idx()), 
  mParticle3Idx(p3Idx), 
  mDcaDaughters12(closePair->dcaDaughters()), 
  mDcaDaughters23(std::numeric_limits<float>::max()), 
  mDcaDaughters31(std::numeric_limits<float>::max()),
  mDV0Max(std::numeric_limits<float>::quiet_NaN())
{
  if(!closePair)
  {
    cerr << "StHFTriplet::StHFTriplet(closePair* ..): closepair not found" << endl;
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    mParticle3Idx = std::numeric_limits<unsigned short>::max();
    return;
  }
  if (!particle3 || 
      (closePair->particle1Idx() == particle3->id() || closePair->particle2Idx() == particle3->id())) {
    mParticle1Idx = std::numeric_limits<unsigned short>::max();
    mParticle2Idx = std::numeric_limits<unsigned short>::max();
    mParticle3Idx = std::numeric_limits<unsigned short>::max();
    return;
  }

  StPhysicalHelixD p3Helix = particle3->helix(bField);

  calculateTopology(closePair, p3Helix, p3MassHypo, particle3->charge(), p3Idx, vtx, bField);
}


// _________________________________________________________
void StHFTriplet::calculateTopology(StHFClosePair * closePair, StPhysicalHelixD & p3Helix, 
				    float p3MassHypo, float p3Charge,
				    unsigned short p3Idx,
				    StThreeVectorF const & vtx, float bField)
{

  mParticle1Dca = closePair->particle1Dca();
  mParticle2Dca = closePair->particle2Dca();
  mParticle1Idx = closePair->particle1Idx();
  mParticle2Idx = closePair->particle2Idx();
  mParticle3Idx = p3Idx;
  mDcaDaughters12 = closePair->dcaDaughters();

  p3Helix.moveOrigin(p3Helix.pathLength(vtx));

  StThreeVectorF const p3Mom = p3Helix.momentum(bField * kilogauss);

  StPhysicalHelixD * p1StraightLine = closePair->p1StraightLine();
  StPhysicalHelixD * p2StraightLine = closePair->p2StraightLine();

  if (!p1StraightLine || !p2StraightLine)
  {
    cerr << "StHFTriplet::StHFTriplet(closePair* ...): Straightlines not found in StHFClosePair" << endl;
    throw;
  }

  StPhysicalHelixD const p3StraightLine(p3Mom, p3Helix.origin(), 0, p3Charge);

  StThreeVectorF p1AtDcaToP2 = closePair->p1AtDcaToP2();
  StThreeVectorF p2AtDcaToP1 = closePair->p2AtDcaToP1();

  std::pair<double, double> const ss23 = p2StraightLine->pathLengths(p3StraightLine);
  StThreeVectorF const p2AtDcaToP3 = p2StraightLine->at(ss23.first);
  StThreeVectorF const p3AtDcaToP2 = p3StraightLine.at(ss23.second);
  
  std::pair<double, double> const ss13 = p1StraightLine->pathLengths(p3StraightLine);
  StThreeVectorF const p1AtDcaToP3 = p1StraightLine->at(ss13.first);
  StThreeVectorF const p3AtDcaToP1 = p3StraightLine.at(ss13.second);

  // -- calculate DCA of particle2 to particl3 at their DCA
  mDcaDaughters23 = (p2AtDcaToP3 - p3AtDcaToP2).mag();
  
  // -- calculate DCA of particle3 to particle1 at their DCA
  mDcaDaughters31 = (p3AtDcaToP1 - p1AtDcaToP3).mag();

  // -- calculate vertices between particles
  StThreeVectorF decayVertex12 = (p1AtDcaToP2 + p2AtDcaToP1) / 2.;
  StThreeVectorF decayVertex23 = (p2AtDcaToP3 + p3AtDcaToP2) / 2.;
  StThreeVectorF decayVertex31 = (p3AtDcaToP1 + p1AtDcaToP3) / 2.;
  
  // -- calculate decay vertex (secondary)
  StThreeVectorF mDecayVertex = ( decayVertex12 + decayVertex23 + decayVertex31) / 3.0;

  // -- calculate maximum distance between vertices of pairs
  const float vDist1223 = (decayVertex12 - decayVertex23).mag();
  const float vDist2331 = (decayVertex23 - decayVertex31).mag();
  const float vDist3112 = (decayVertex31 - decayVertex12).mag();
  mDV0Max = vDist1223 > vDist2331 ? vDist1223 : vDist2331;
  mDV0Max = mDV0Max > vDist3112 ? mDV0Max : vDist3112;
  
  // -- constructing mother daughter four momentum. Need helix (not straight line) for each daughter
  if (!closePair->p1Helix() || !closePair->p2Helix())
  {
    cerr << "StHFTriplet::StHFTriplet(closePair* ...): Helices not found in StHFClosePair" << endl;
    throw;
  }

  double const p1AtV0 = closePair->p1Helix()->pathLength( mDecayVertex );
  StThreeVectorF const p1MomAtDca = closePair->p1Helix()->momentumAt(p1AtV0 ,  bField * kilogauss);

  double const p2AtV0 = closePair->p2Helix()->pathLength( mDecayVertex );
  StThreeVectorF const p2MomAtDca = closePair->p2Helix()->momentumAt(p2AtV0 ,  bField * kilogauss);
  
  double const p3AtV0 = p3Helix.pathLength( mDecayVertex );
  StThreeVectorF const p3MomAtDca = p3Helix.momentumAt(p3AtV0 ,  bField * kilogauss);
  
  StLorentzVectorF const p1FourMom(p1MomAtDca, p1MomAtDca.massHypothesis(closePair->p1massHypothesis()));
  StLorentzVectorF const p2FourMom(p2MomAtDca, p2MomAtDca.massHypothesis(closePair->p2massHypothesis()));
  StLorentzVectorF const p3FourMom(p3MomAtDca, p3MomAtDca.massHypothesis(p3MassHypo));
  
  mLorentzVector = p1FourMom + p2FourMom + p3FourMom;
   
  // -- calculate pointing angle and decay length
  StThreeVectorF const vtxToV0 = mDecayVertex - vtx;
  mPointingAngle = vtxToV0.angle(mLorentzVector.vect());
  mDecayLength = vtxToV0.mag();

}

#ifndef StHFTriplet_hh
#define StHFTriplet_hh

/* **************************************************
 *  Generic class calculating and storing primary triplets in HF analysis
 *  Allows to combine:
 *  - three particles, using
 *      StHFTriplet(StPicoTrack const * particle1, StPicoTrack const * particle2, 
 *                  StPicoTrack const * particle3, ...
 *
 * **************************************************
 *
 *  Initial Authors:                                                                                                                       
 *            Xin Dong        (xdong@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)
 *          **Michael Lomnitz (mrlomnitz@lbl.gov)
 *
 *  ** Code Maintainer    
 *
 * **************************************************
 */

#include "TObject.h"
#include "TClonesArray.h"

#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"

class StPicoTrack;
class StPicoEvent;
class StHFClosePair;

class StHFTriplet : public TObject
{
 public:
  StHFTriplet();
  explicit StHFTriplet(StHFTriplet const *);
  StHFTriplet(StPicoTrack const * particle1, StPicoTrack const * particle2, StPicoTrack const * particle3, 
	     float p1MassHypo, float p2MassHypo, float p3MassHypo,
	     unsigned short p1Idx, unsigned short p2Idx, unsigned short p3Idx,
	     StThreeVectorF const & vtx, float bField);
  StHFTriplet(StHFClosePair * pair, StPicoTrack const * particle3, 
	      float p3MassHypo,
	      unsigned short p3Idx,
	      StThreeVectorF const & vtx, float bField);
  ~StHFTriplet() {;}

  StLorentzVectorF const & lorentzVector() const;
  StThreeVectorF const & decayVertex() const;
  inline float m()    const;
  inline float pt()   const;
  inline float eta()  const;
  inline float phi()  const;
  inline float pointingAngle() const;
  inline float decayLength() const;
  inline float particle1Dca() const;
  inline float particle2Dca() const;
  inline float particle3Dca() const;
  inline unsigned short particle1Idx() const;
  inline unsigned short particle2Idx() const;
  inline unsigned short particle3Idx() const;
  inline float dcaDaughters12() const;
  inline float dcaDaughters23() const;
  inline float dcaDaughters31() const;
  inline float v0x() const;
  inline float v0y() const;
  inline float v0z() const;
  inline float px() const;
  inline float py() const;
  inline float pz() const;
  inline float DcaToPrimaryVertex() const;
  inline float dV0Max() const;

 protected:
  void calculateTopology(StHFClosePair * pair, StPhysicalHelixD  & p3Helix, 
			 float p3MassHypo, float p3Charge,
			 unsigned short p3Idx,
			 StThreeVectorF const & vtx, float bField);
 private:
  StHFTriplet(StHFTriplet const &);
  StHFTriplet& operator=(StHFTriplet const &);
  StLorentzVectorF mLorentzVector; 
  StThreeVectorF   mDecayVertex; 

  float mPointingAngle;
  float mDecayLength;

  float mParticle1Dca;
  float mParticle2Dca;
  float mParticle3Dca;

  unsigned short  mParticle1Idx; // index of track in StPicoDstEvent
  unsigned short  mParticle2Idx;
  unsigned short  mParticle3Idx;

  float mDcaDaughters12;
  float mDcaDaughters23;
  float mDcaDaughters31;

  float mDV0Max;

  ClassDef(StHFTriplet,2)
};
inline StLorentzVectorF const & StHFTriplet::lorentzVector() const { return mLorentzVector;}
inline float StHFTriplet::m()    const { return mLorentzVector.m();}
inline float StHFTriplet::pt()   const { return mLorentzVector.perp();}
inline float StHFTriplet::eta()  const { return mLorentzVector.pseudoRapidity();}
inline float StHFTriplet::phi()  const { return mLorentzVector.phi();}
inline float StHFTriplet::px()   const { return mLorentzVector.px();}
inline float StHFTriplet::py()   const { return mLorentzVector.py();}
inline float StHFTriplet::pz()   const { return mLorentzVector.pz();}
inline float StHFTriplet::pointingAngle() const { return mPointingAngle;}
inline float StHFTriplet::decayLength()   const { return mDecayLength;}
inline float StHFTriplet::particle1Dca()  const { return mParticle1Dca;}
inline float StHFTriplet::particle2Dca()  const { return mParticle2Dca;}
inline float StHFTriplet::particle3Dca()  const { return mParticle3Dca;}
inline unsigned short StHFTriplet::particle1Idx() const { return mParticle1Idx;}
inline unsigned short StHFTriplet::particle2Idx() const { return mParticle2Idx;}
inline unsigned short StHFTriplet::particle3Idx() const { return mParticle3Idx;}
inline float StHFTriplet::dcaDaughters12() const { return mDcaDaughters12;}
inline float StHFTriplet::dcaDaughters23() const { return mDcaDaughters23;}
inline float StHFTriplet::dcaDaughters31() const { return mDcaDaughters31;}
inline StThreeVectorF const & StHFTriplet::decayVertex() const { return mDecayVertex;}
inline float StHFTriplet::v0x() const { return mDecayVertex.x();}
inline float StHFTriplet::v0y() const { return mDecayVertex.y();}
inline float StHFTriplet::v0z() const { return mDecayVertex.z();}
inline float StHFTriplet::DcaToPrimaryVertex() const {return mDecayLength*sin(mPointingAngle);}
inline float StHFTriplet::dV0Max() const {return mDV0Max;}
#endif

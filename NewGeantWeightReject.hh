bool NewGeantWeightReject ( const double pt, const double weight, const int type ){
  // Reject jets from pT hard bins with very few entries and consequently a high weight
  // Couldn't find a better way than compiling a lookup table by hand.
  // IMPORTANT! Jet pT depends on R (and other jet finding parameters,
  // but they should all be unchanged for GEANT MC data)
  // --> You should always supply all info and not be lulled by default values
  // If R changes, you need to get a fresh lookup table!

  // type:
  // MC: 0 for trigger, 1 for recoil, 2 for inclusive.
  // RECO: 10 for trigger, 11 for recoil, 12 for inclusive.
  // NOTE: Determined with nominal cuts. May need to adjust
  
  // cout << "R = " << R << endl;
  // cout << "type = " << type << endl;
  // cout << "pt = " << pt << endl;
  // cout << "weight = " << weight << endl;
    
  // Note: CINT can be dumb with switch statements, so I leave in break's.

  if ( 1 ) {
    // Always accept pT>40. You shouldn't be using it, but if you do, use whatever you got.
    if ( pt>=40 ) return false;

    // Always accept pT<5. They're all okay
    if ( pt<5 ) return false;
    
    switch (type){
    case 0:  // trigger jets
      if ( ( 5 < pt && pt < 10) )               return false; // 100s of cases
      if ( (10 < pt && pt < 15) && weight>400 ) return true; // 22 cases (14 in |eta|<1)
      if ( (15 < pt && pt < 20) && weight>50  ) return true; // 21 cases (19 in |eta|<1)
      if ( (20 < pt && pt < 25) && weight>13  ) return true; // 10 cases (10 in |eta|<1)
      if ( (25 < pt && pt < 30) && weight>5   ) return true; // 6 cases (6 in |eta|<1)
      // if ( (30 < pt && pt < 35) && weight>0.9 ) return true; // 46 cases (43 in |eta|<1)
      if ( (35 < pt && pt < 40) && weight>0.9 ) return true; // 5 cases (5 in |eta|<1)
      return false;
      break;
    case 1:  // recoil jets
      if ( ( 5 < pt && pt < 10) && weight>1800    ) return true; // lots of cases (13 in |eta|<1)
      if ( (10 < pt && pt < 15) && weight>50    ) return true; // 48 cases (14 in |eta|<1)
      if ( (15 < pt && pt < 20) && weight>5     ) return true; // lots of cases (19 in |eta|<1)
      if ( (20 < pt && pt < 25) && weight>5     ) return true; // 3 cases (1 in |eta|<1)
      if ( (25 < pt && pt < 30) && weight>0.9   ) return true; // 4 cases (4 in |eta|<1)
      // if ( (30 < pt && pt <35) && weight>0.9 ) return true; // 100s of cases (100s in |eta|<1)
      if ( (35 < pt && pt <40) && weight>0.025 ) return true; // 10 cases (10 in |eta|<1)
      return false;
      break;
    case 2:  // inclusive jets
      if ( ( 5 < pt && pt < 10) )                return false; // ~1500 cases, good enough
      if ( (10 < pt && pt < 15) && weight>300 )  return true; // ~100 cases
      if ( (15 < pt && pt < 20) && weight>50  )  return true; // 17 cases
      if ( (20 < pt && pt < 25) && weight>5  )   return true; // ~100 cases
      if ( (25 < pt && pt < 30) && weight>5  )   return true; // 6 cases. 
      if ( (30 < pt && pt < 35) && weight>0.9 ) return true; // 40 cases
      if ( (35 < pt && pt < 40) && weight>0.9 ) return true; // 4 cases
      return false;
      break;
    case 10:  // trigger jets - RECO
      if ( pt < 10 )               return false; // should never happen anyway
      if ( (10 < pt && pt < 15) && weight>50 ) return true; // 17 cases
      if ( (15 < pt && pt < 20) && weight>13 ) return true; // 9 cases
      if ( (20 < pt && pt < 25) && weight>5  ) return true; // 11 cases
      // if ( (25 < pt && pt < 30) && weight>5  ) return true; // 70 cases
      if ( (30 < pt && pt < 35) && weight>0.9 ) return true; // 7 cases
      // if ( (35 < pt && pt < 40) && weight>0.9 ) return true; // 181 cases
      return false;
      break;
    case 11:  // recoil jets - RECO
      if ( pt < 10 )               return false; // should never happen anyway
      if ( (10 < pt && pt < 15) && weight>13    ) return true; // 1 case
      if ( (15 < pt && pt < 20) && weight>5     ) return true; // 2 cases
      // if ( (20 < pt && pt < 25) && weight>5     ) return true; // 82 cases
      if ( (25 < pt && pt < 30) && weight>0.9   ) return true; // 3 cases 
      // if ( (30 < pt && pt <35) && weight>0.9 ) return true; // 100s of cases
      if ( (35 < pt && pt <40) && weight>0.025 ) return true; // 10 cases - borderline
      return false;
      break;
    case 12:  // inclusive jets
      if ( ( 5 < pt && pt < 10) && weight>200)   return true; // ~100 cases
      if ( (10 < pt && pt < 15) && weight>40 )   return true; // ~50 cases
      if ( (15 < pt && pt < 20) && weight>10  )  return true; // ~12 cases
      if ( (20 < pt && pt < 25) && weight>5  )   return true; // ~40 cases
      if ( (25 < pt && pt < 30) && weight>5   )  return true; // 1 case. may want to go down to 0.9 and exclude 100 more
      if ( (30 < pt && pt < 35) && weight>0.9 )  return true; // 9 cases
      if ( (35 < pt && pt < 40) && weight>0.9 )  return true; // 0 cases actually
      return false;
      break;
    default:
      cerr << "Unknown type " << type << "! Bailing out." << endl;
      return true;
      break;
    }

  }
  
  // If no acceptable R, reject
  return true;
}
LorentzVectorBase.coordinate_names(::LorentzVectorBase.PtEtaPhiM) = (:pt, :eta, :phi, :mass)

LorentzVectorBase.t(::LorentzVectorBase.PtEtaPhiM,lv) = LorentzVectorBase.energy(LorentzVectorBase.PtEtaPhiM(),lv)
LorentzVectorBase.x(::LorentzVectorBase.PtEtaPhiM,lv) = LorentzVectorBase.px(lv)
LorentzVectorBase.y(::LorentzVectorBase.PtEtaPhiM,lv) = LorentzVectorBase.py(lv)
LorentzVectorBase.z(::LorentzVectorBase.PtEtaPhiM,lv) = LorentzVectorBase.pz(lv)


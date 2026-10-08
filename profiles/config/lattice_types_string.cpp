//
// Created by genshen on 2019/9/19.
//

#include "lattice_types_string.h"
#include <stdexcept>

std::string lat::LatTypesString(LatticeTypes::lat_type lat_type) {
  switch (lat_type) {
    case LatticeTypes::V:
      return "V";
    case LatticeTypes::Mo:
      return "Mo";
    case LatticeTypes::Re:
      return "Re";
    case LatticeTypes::Mn:
      return "Mn";
    case LatticeTypes::Ni:
      return "Ni";
    case LatticeTypes::Si:
      return "Si";
    case LatticeTypes::MoMo:
      return "MoMo";
    case LatticeTypes::MoRe:
      return "MoRe";
    case LatticeTypes::MoMn:
      return "MoMn";
    case LatticeTypes::MoNi:
      return "MoNi";
    case LatticeTypes::MoSi:
      return "MoSi";
    case LatticeTypes::ReRe:
      return "ReRe";
    case LatticeTypes::ReMn:
      return "ReMn";
    case LatticeTypes::ReNi:
      return "ReNi";
    case LatticeTypes::ReSi:
      return "ReSi";
    case LatticeTypes::MnMn:
      return "MnMn";
    case LatticeTypes::MnNi:
      return "MnNi";
    case LatticeTypes::MnSi:
      return "MnSi";
    case LatticeTypes::NiNi:
      return "NiNi";
    case LatticeTypes::NiSi:
      return "NiSi";
    case LatticeTypes::SiSi:
      return "SiSi";
  }
}

LatticeTypes::lat_type lat::LatTypes(const std::string &lat_type) {
  if (lat_type == "V") {
    return LatticeTypes::V;
  } else if (lat_type == "Mo") {
    return LatticeTypes::Mo;
  } else if (lat_type == "Re") {
    return LatticeTypes::Re;
  } else if (lat_type == "Mn") {
    return LatticeTypes::Mn;
  } else if (lat_type == "Ni") {
    return LatticeTypes::Ni;
  } else if (lat_type == "Si") {
    return LatticeTypes::Si;
  } else if (lat_type == "MoMo") {
    return LatticeTypes::MoMo;
  } else if (lat_type == "MoRe") {
    return LatticeTypes::MoRe;
  } else if (lat_type == "MoMn") {
    return LatticeTypes::MoMn;
  } else if (lat_type == "MoNi") {
    return LatticeTypes::MoNi;
  } else if (lat_type == "MoSi") {
    return LatticeTypes::MoSi;
  } else if (lat_type == "ReRe") {
    return LatticeTypes::ReRe;
  } else if (lat_type == "ReMn") {
    return LatticeTypes::ReMn;
  } else if (lat_type == "ReNi") {
    return LatticeTypes::ReNi;
  } else if (lat_type == "ReSi") {
    return LatticeTypes::ReSi;
  } else if (lat_type == "MnMn") {
    return LatticeTypes::MnMn;
  } else if (lat_type == "MnNi") {
    return LatticeTypes::MnNi;
  } else if (lat_type == "MnSi") {
    return LatticeTypes::MnSi;
  } else if (lat_type == "NiNi") {
    return LatticeTypes::NiNi;
  } else if (lat_type == "NiSi") {
    return LatticeTypes::NiSi;
  } else if (lat_type == "SiSi") {
    return LatticeTypes::SiSi;
  } else {
    throw std::invalid_argument("wrong lattice type");
  }
}

#pragma once

#include "../placement/placement_analyzer.hpp"

class analyzer_params
{
public:
  bool only_critical = true;
  bool preserve_depth = true;
  bool use_fanouts = true;
  bool preserve_interval = false;
  mockturtle::placement_params placement_ps = {};
};
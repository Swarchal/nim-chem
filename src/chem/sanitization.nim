## The standard clean-up run after parsing.

import ./[mol, valence, aromaticity, stereo]

proc sanitize*(m: var Mol) =
  ## Kekulizes, checks valences on the Kekulé structure, perceives
  ## aromaticity, then removes stereo that is not real (`assignStereo`).
  ## Raises `KekulizeError` or `ValenceError`.
  kekulize(m)
  checkValences(m)
  setAromaticity(m)
  assignStereo(m)

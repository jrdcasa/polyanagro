function output = isOctave
  output = (exist ('OCTAVE_VERSION', 'builtin') > 0);
end
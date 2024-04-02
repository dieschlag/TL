% Stéphane Rossignol - 2021 pour CentraleSupélec

function [] = flush(num)

aa=version('-release'); %%% si c'est de longueur nulle, c'est 'octave' qui est utilisé (est-ce sûr ???)

if length(aa)==0  %%% ça veut dire qu'on utilise 'octave' (bien sûr, ce n'est pas bien fait)
  fflush(num);    %%% => 'fflush' n'existe pas sous 'matlab' :-(
end;


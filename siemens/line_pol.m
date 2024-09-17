%-------------------------------------------------------------------------------
% find out the polarity of the line
%-------------------------------------------------------------------------------
function y=line_pol(iLineOfContr,iContr,para)
  b_epi_positive=para.b_epi_positive;
  b_epi_pol=para.b_epi_pol;
  if b_epi_positive
      % all lines are either positive or negative
      if b_epi_pol(iContr) 
          y=1;
      else
          y=0;
      end
  else
      % decide the polarity of the line based on the epi polarity and the index of the line
      if (b_epi_pol(iContr) && (mod(iLineOfContr-1,2) == 0)) || ...
              (~b_epi_pol(iContr) && (mod(iLineOfContr-1, 2) == 1))
          y=1;
      else
          y=0;
      end
  end
end

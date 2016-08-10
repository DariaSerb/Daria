function out = Validation(P)
      
  if (P.radiusout < 0.5 * P.lx)&&(P.radiusout < 0.5 * P.ly)   
  out = true;
  disp('correct initial parameters');
  else
  out = false;
  disp('uncorrect initial parameters'); 
  end
end

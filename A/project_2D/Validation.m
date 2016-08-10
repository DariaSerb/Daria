function out = Validation(P)
      
  if (P.radiusout < 0.5*P.lx)&&(P.radiusout < 0.5*P.ly)   
  out = true;
  else
  out = false;    
  end
end

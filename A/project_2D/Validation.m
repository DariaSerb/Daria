function out = Validation(P)
      
  if (P.radiusout < 0.5 * P.lx)&&(P.radiusout < 0.5 * P.ly)   
  out = true;
  disp('correct input data');
  else
  out = false;
  disp('incorrect input data: to change the size of mesh or delta'); 
  end
end

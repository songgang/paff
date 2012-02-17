function varScope3 
    a = @(x) nestfun1(x, 3);
    a(1)
    
   function z = nestfun1(x, y)
      z = x+y
   
 

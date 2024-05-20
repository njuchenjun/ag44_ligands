  Program FBC

  Use ctrl

  Implicit None 

  Call INIT_PES() !! jyh 22/9/6

  Call Initialize()
  Call simulation()
  
  Call DEALLOCATE_PES !! jyh 22/9/6

  End Program

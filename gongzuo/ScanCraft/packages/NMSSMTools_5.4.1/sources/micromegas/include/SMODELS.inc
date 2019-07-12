{  
   char smodels_database[20];
   char smodels_version[20];
   printf("\nSMODELS:\n");
   smodels(5/*nf*/, 0.1/*fb*/, "smodels.in",0);
   system("make -C ../Packages/smodels 1>/dev/null"); 
   system("../Packages/smodels/runTools.py -particle ./ xseccomputer -p -N -O -f smodels.in");
   system("../Packages/smodels/runSModelS.py -f smodels.in -s smodels.res -particles ./  > smodels.out "); 
   slhaRead("smodels.res", 1);
   slhaSTRFormat("SModelS_Settings"," 0 %s",smodels_version);
   slhaSTRFormat("SModelS_Settings"," 1 %s",smodels_database);
   printf("%s with database %s \n",smodels_version, smodels_database);
   result=slhaVal("SModelS_Exclusion",0.,2,0,0); 
   if(result>=0)
   { 
     slhaSTRFormat("SModelS_Exclusion"," 1 0  %s",topology);
     Rvalue=slhaVal("SModelS_Exclusion",0.,2,1,1);
     slhaSTRFormat("SModelS_Exclusion"," 1 3 %s",analysis);
     printf("highest R=%.2E from %s, topology %s\n", Rvalue,analysis,topology);
   }
   switch(result)
   { case -1: printf("no channels for testing.\n");break;
     case  0: printf("not excluded.\n");break; 
     case  1: printf("excluded.\n ");break;
   }  
}
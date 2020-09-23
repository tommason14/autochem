from autochem import Settings              
                                           
sett=Settings()                            
sett.input.contrl.runtyp='energy'          
sett.input.basis.gbasis='cct'              
sett.input.mp2.scsopo=1.64                 
sett.input.fmo.nbody=3                     
# trimer cutoffs giving < 2 kJ/mol error cf. full system
sett.input.fmo.ritrim='1.77,1.77,1.77,1.77'

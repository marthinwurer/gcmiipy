



"""
Take the maps from 611 in hansen 1983 and use them as a baseline
"""

land_cover_map = \
"""
       01475667840 120     02       
0 01  4535623AXX4   10 122478843420 
628A99898737187330 3896799AXXXXXXXA9
005349AAA337211  02759AXXXXXXAXX6452
     2AXA99A61   169AXXXXXXXXAXX6 2 
      9XX9951    37586769XXXXX973   
      5XXX8      37623A9XXXXXX532   
      07922      8XXAA988AXXXX6     
  0    07411    2XXXXX7A3394A41     
         3342   1AXAXX94 04 4310    
          3X93   237XXA2  0 4241    
          6XXA6    4X95     2422510 
          2XXX8    2XA51     01431  
0          6XX5    2XX34     05X93 0
           6X6     09801     1XXA8  
1          891      31        4175 0
1          81                    1 2
           60            0          
           0                        
           12     00012543356545531 
000134544769810147AXXXXXXXXXXXXXXX93
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"""

land_cover_key = {
    " ": 0,
    "0": 0.025,
    "1": 0.1,
    "2": 0.2,
    "3": 0.3,
    "4": 0.4,
    "5": 0.5,
    "6": 0.6,
    "7": 0.7,
    "8": 0.8,
    "9": 0.9,
    "A": 0.975,
    "X": 1,
}

topography_map = \
"""
       0012237AA50 000     00       
0 00  0101102CPP6   00 000011110000 
30278523110208H610 23201211333235544
001369753002111  0020122322346654321
     1A8421320   0134221234BC9742 1 
      AI72210    133327319HDHD642   
      5G722      262037CAK+++A111   
      07700      346337575BON62     
  0    08100    024565661142611     
         2120   033455A3 02 1100    
          3430   0045781  0 3220    
          52221    15A5     2120610 
          4B243    1BA20     00101  
0          J463    2B712     02421 0
           C23     09800     04323  
0          820      21        1013 0
0          50                    0 1
           2             0          
           0
           01     00002761489A88720 
00002675467960004BLTTSUUNITXWUTRPME2
543347BEIKHA7778FMRVXYZZY+++YWUSPNG8
RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"""

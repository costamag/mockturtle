module top( x0 , x1 , x2 , x3 , x4 , y0 , y1 );
  input x0 , x1 , x2 , x3 , x4 ;
  output y0 , y1 ;
  wire n0 , n2 , n7 , n8 , n10 ;
  _const0_             g0( .z (n0) );
  _const1_             g1( .z (n2) );
  AND2x2_ASAP7_75t_R   g2( .A (x1), .B (x3), .Y (n7) );
  INVx1_ASAP7_75t_R    g3( .A (n7), .Y (n8) );
  AO22x1_ASAP7_75t_R   g4( .A1 (n8), .A2 (x2), .B1 (x1), .B2 (x0), .Y (y0) );
  NAND2x1_ASAP7_75t_R  g5( .A (x3), .B (x1), .Y (n10) );
  OA21x2_ASAP7_75t_R   g6( .A1 (x4), .A2 (x2), .B (n10), .Y (y1) );
endmodule

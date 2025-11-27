module top( x0 , x1 , x2 , x3 , x4 , y0 , y1 );
  input x0 , x1 , x2 , x3 , x4 ;
  output y0 , y1 ;
  wire n0 , n2 , n7 , n8 , n9 , n11 ;
  _const0_             g0( .z (n0) );
  _const1_             g1( .z (n2) );
  INVx1_ASAP7_75t_R    g2( .A (x2), .Y (n7) );
  AND2x2_ASAP7_75t_R   g3( .A (x1), .B (x3), .Y (n8) );
  NAND2x1_ASAP7_75t_R  g4( .A (x0), .B (x1), .Y (n9) );
  OAI21x1_ASAP7_75t_R  g5( .A1 (n7), .A2 (n8), .B (n9), .Y (y0) );
  INVx1_ASAP7_75t_R    g6( .A (x4), .Y (n11) );
  AOI22x1_ASAP7_75t_R  g7( .A1 (x1), .A2 (x3), .B1 (n11), .B2 (n7), .Y (y1) );
endmodule

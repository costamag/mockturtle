module top( x0 , x1 , x2 , x3 , x4 , y0 , y1 );
  input x0 , x1 , x2 , x3 , x4 ;
  output y0 , y1 ;
  wire n0 , n2 , n7 ;
  _const0_             g0( .z (n0) );
  _const1_             g1( .z (n2) );
  NAND2x1_ASAP7_75t_R  g2( .A (x1), .B (x3), .Y (n7) );
  AO22x1_ASAP7_75t_R   g3( .A1 (x0), .A2 (x1), .B1 (n7), .B2 (x2), .Y (y0) );
  OA21x2_ASAP7_75t_R   g4( .A1 (x4), .A2 (x2), .B (n7), .Y (y1) );
endmodule

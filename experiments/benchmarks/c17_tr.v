// Benchmark "c17" written by ABC on Sun Dec 22 23:09:28 2024

module c17 ( 
    G1, G3, G2, G4, G5,
    G16, G17  );
  input  G1, G3, G2, G4, G5;
  output G16, G17;
  wire new_n8, new_n9, new_n10, new_n12;
  INVx1_ASAP7_75t_R    g0(.A(G2), .Y(new_n8));
  NAND2x1_ASAP7_75t_R  g1(.A(G1), .B(G3), .Y(new_n9));
  AND2x2_ASAP7_75t_R   g2(.A(G3), .B(G4), .Y(new_n10));
  OAI21x1_ASAP7_75t_R  g3(.A1(new_n8), .A2(new_n10), .B(new_n9), .Y(G16));
  INVx1_ASAP7_75t_R    g4(.A(G5), .Y(new_n12));
  AOI22x1_ASAP7_75t_R  g5(.A1(G3), .A2(G4), .B1(new_n8), .B2(new_n12), .Y(G17));
endmodule



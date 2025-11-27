// Benchmark "top" written by ABC on Sun Dec 22 23:24:08 2024

module top ( 
    x0, x1, x2, x3, x4,
    y0, y1  );
  input  x0, x1, x2, x3, x4;
  output y0, y1;
  wire new_n8;
  assign new_n8 = ~x1 | ~x3;
  assign y0 = (new_n8 & x2) | (x0 & x1);
  assign y1 = new_n8 & (x2 | x4);
endmodule



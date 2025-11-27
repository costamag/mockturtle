// Benchmark "top" written by ABC on Sun Dec 22 23:30:00 2024

module top ( 
    x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15,
    x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29,
    x30, x31, x32, x33, x34, x35,
    y0, y1, y2, y3, y4, y5, y6  );
  input  x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14,
    x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28,
    x29, x30, x31, x32, x33, x34, x35;
  output y0, y1, y2, y3, y4, y5, y6;
  wire new_n44, new_n45, new_n46, new_n47, new_n48, new_n49, new_n50,
    new_n51, new_n52, new_n53, new_n54, new_n55, new_n56, new_n57, new_n58,
    new_n60, new_n61, new_n62, new_n63, new_n64, new_n65, new_n66, new_n67,
    new_n68, new_n69, new_n70, new_n71, new_n72, new_n73, new_n74, new_n75,
    new_n76, new_n77, new_n79, new_n80, new_n81, new_n82, new_n83, new_n84,
    new_n85, new_n86, new_n87, new_n88, new_n89, new_n90, new_n91, new_n92,
    new_n93, new_n94, new_n95, new_n96, new_n97, new_n98, new_n99,
    new_n100, new_n102, new_n103, new_n104, new_n105, new_n106, new_n107,
    new_n108, new_n109, new_n110, new_n111, new_n112, new_n113, new_n114,
    new_n117, new_n118, new_n119, new_n121, new_n122, new_n123;
  assign new_n44 = ~x2;
  assign new_n45 = ~x16;
  assign new_n46 = ~x5;
  assign new_n47 = ~x4 & ~new_n46;
  assign new_n48 = ~new_n47 & (~new_n45 | ~x17) & (~new_n44 | ~x3);
  assign new_n49 = ~x11;
  assign new_n50 = ~x13;
  assign new_n51 = (new_n49 | x10) & (new_n50 | x12);
  assign new_n52 = ~x7;
  assign new_n53 = ~x15;
  assign new_n54 = ~x0;
  assign new_n55 = ~x8;
  assign new_n56 = (~new_n55 | ~x9) & (~new_n54 | ~x1);
  assign new_n57 = new_n56 & (x6 | new_n52) & (x14 | new_n53);
  assign new_n58 = new_n57 & new_n48 & new_n51;
  assign y0 = ~new_n58;
  assign new_n60 = ~x18;
  assign new_n61 = ~x19;
  assign new_n62 = ~new_n46 & ~x20 & (~y0 | ~x4);
  assign new_n63 = (new_n58 | new_n44) & x3;
  assign new_n64 = (new_n58 | new_n54) & x1;
  assign new_n65 = ~new_n62 & (~new_n61 | ~new_n63) & (~new_n60 | ~new_n64);
  assign new_n66 = ~x9;
  assign new_n67 = new_n66 | x22 | (x8 & y0);
  assign new_n68 = new_n53 | x25 | (x14 & y0);
  assign new_n69 = new_n52 | x21 | (x6 & y0);
  assign new_n70 = new_n50 | x24 | (x12 & y0);
  assign new_n71 = ~x26;
  assign new_n72 = (new_n58 | new_n45) & x17;
  assign new_n73 = ~new_n71 | ~new_n72;
  assign new_n74 = new_n49 | (y0 & x10);
  assign new_n75 = (x23 | new_n74) & new_n73;
  assign new_n76 = new_n70 & new_n69 & new_n68 & new_n67 & new_n75;
  assign new_n77 = new_n76 & new_n65;
  assign y1 = ~new_n77;
  assign new_n79 = new_n46 | (x20 & y1) | (x4 & y0);
  assign new_n80 = (~new_n61 & ~new_n77) | ~new_n63;
  assign new_n81 = (~new_n60 & ~new_n77) | ~new_n64;
  assign new_n82 = ~x27 & ~new_n81;
  assign new_n83 = ~new_n82;
  assign new_n84 = new_n83 & (x28 | new_n80) & (x29 | new_n79);
  assign new_n85 = ~new_n53 & (~x14 | ~y0);
  assign new_n86 = ~new_n85;
  assign new_n87 = x34 | new_n86 | (new_n68 & y1);
  assign new_n88 = x33 | new_n50;
  assign new_n89 = new_n88 | (new_n70 & y1) | (x12 & y0);
  assign new_n90 = new_n66 | x31 | (x8 & y0);
  assign new_n91 = new_n90 | (y1 & new_n67);
  assign new_n92 = ~new_n52 & (~x6 | ~y0);
  assign new_n93 = ~new_n92;
  assign new_n94 = x30 | new_n93 | (new_n69 & y1);
  assign new_n95 = ~x35;
  assign new_n96 = ~x32 & ~new_n74 & (~y1 | ~x23);
  assign new_n97 = (new_n77 | new_n71) & new_n72;
  assign new_n98 = ~new_n96 & (~new_n95 | ~new_n97);
  assign new_n99 = new_n94 & new_n91 & new_n89 & new_n87 & new_n98;
  assign new_n100 = new_n99 & new_n84;
  assign y2 = ~new_n100;
  assign new_n102 = ~x31 | ~y2;
  assign new_n103 = ~new_n66 & (~x22 | ~y1) & (~x8 | ~y0);
  assign new_n104 = ~new_n74 & (~x32 | ~y2) & (~x23 | ~y1);
  assign new_n105 = ~new_n86 & (~x34 | ~y2) & (~x25 | ~y1);
  assign new_n106 = new_n50 | (x24 & y1) | (x12 & y0);
  assign new_n107 = ~new_n106 & (~x33 | ~y2);
  assign new_n108 = new_n107 | new_n105;
  assign new_n109 = ~new_n108 & ~new_n104 & (~new_n102 | ~new_n103);
  assign new_n110 = new_n80 | (y2 & x28);
  assign new_n111 = new_n93 | (x30 & y2) | (x21 & y1);
  assign new_n112 = new_n79 | (y2 & x29);
  assign new_n113 = new_n112 & new_n110 & new_n83 & new_n81 & new_n111;
  assign new_n114 = new_n97 & new_n99 & new_n84;
  assign y3 = ~new_n114 & (~new_n113 | ~new_n109) & (~new_n95 | ~new_n97);
  assign y4 = ~new_n109;
  assign new_n117 = ~new_n104 & ~new_n112 & (~new_n102 | ~new_n103);
  assign new_n118 = new_n107 | new_n111;
  assign new_n119 = ~new_n118 & ~new_n104 & (~new_n102 | ~new_n103);
  assign y5 = new_n119 | new_n108 | new_n117;
  assign new_n121 = ~new_n107;
  assign new_n122 = ~new_n110;
  assign new_n123 = new_n104 | (new_n122 & new_n112);
  assign y6 = new_n105 | new_n119 | (new_n121 & new_n123);
endmodule



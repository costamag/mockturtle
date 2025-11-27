// Benchmark "top" written by ABC on Sun Dec 22 23:47:43 2024

module top ( 
    x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15,
    x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29,
    x30, x31, x32, x33, x34, x35,
    y0, y1, y2, y3, y4, y5, y6  );
  input  x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14,
    x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28,
    x29, x30, x31, x32, x33, x34, x35;
  output y0, y1, y2, y3, y4, y5, y6;
  wire new_n38, new_n39, new_n40, new_n41, new_n42, new_n43, new_n44,
    new_n45, new_n46, new_n47, new_n48, new_n49, new_n50, new_n51, new_n52,
    new_n54, new_n55, new_n56, new_n57, new_n58, new_n59, new_n60, new_n61,
    new_n62, new_n63, new_n64, new_n65, new_n67, new_n68, new_n69, new_n70,
    new_n71, new_n72, new_n73, new_n74, new_n75, new_n76, new_n77, new_n78,
    new_n79, new_n80, new_n81, new_n82, new_n83, new_n84, new_n85, new_n86,
    new_n87, new_n89, new_n90, new_n91, new_n92, new_n93, new_n94, new_n95,
    new_n96, new_n97, new_n98, new_n99, new_n100, new_n101, new_n102,
    new_n103, new_n104, new_n106, new_n107, new_n108, new_n109, new_n110,
    new_n111, new_n112, new_n113, new_n115, new_n116, new_n117, new_n118,
    new_n120, new_n121;
  assign new_n38 = ~x7;
  assign new_n39 = ~x9;
  assign new_n40 = ~x13;
  assign new_n41 = ~x11;
  assign new_n42 = (new_n41 | x10) & (new_n40 | x12);
  assign new_n43 = (new_n39 | x8) & (x6 | new_n38) & new_n42;
  assign new_n44 = ~x17;
  assign new_n45 = x16 | new_n44;
  assign new_n46 = ~x15;
  assign new_n47 = x14 | new_n46;
  assign new_n48 = ~x5;
  assign new_n49 = ~x3;
  assign new_n50 = ~x1;
  assign new_n51 = (x2 | new_n49) & (new_n48 | x4) & (x0 | new_n50);
  assign new_n52 = new_n51 & new_n47 & new_n43 & new_n45;
  assign y0 = ~new_n52;
  assign new_n54 = new_n46 | x25 | (new_n47 & y0);
  assign new_n55 = new_n44 | x26 | (new_n45 & y0);
  assign new_n56 = new_n54 & new_n55;
  assign new_n57 = new_n41 | x23 | (x10 & y0);
  assign new_n58 = new_n40 | x24 | (x12 & y0);
  assign new_n59 = new_n38 | x21 | (x6 & y0);
  assign new_n60 = new_n39 | x22 | (x8 & y0);
  assign new_n61 = new_n60 & new_n59 & new_n58 & new_n56 & new_n57;
  assign new_n62 = new_n48 | x20 | (x4 & y0);
  assign new_n63 = new_n49 | x19 | (x2 & y0);
  assign new_n64 = new_n50 | x18 | (x0 & y0);
  assign new_n65 = new_n64 & new_n63 & new_n61 & new_n62;
  assign y1 = ~new_n65;
  assign new_n67 = new_n41 | x32 | (x10 & y0);
  assign new_n68 = new_n67 | (y1 & new_n57);
  assign new_n69 = new_n40 | x33 | (x12 & y0);
  assign new_n70 = new_n69 | (y1 & new_n58);
  assign new_n71 = new_n38 | x30 | (x6 & y0);
  assign new_n72 = new_n71 | (y1 & new_n59);
  assign new_n73 = new_n39 | x31 | (x8 & y0);
  assign new_n74 = new_n73 | (y1 & new_n60);
  assign new_n75 = new_n74 & new_n72 & new_n68 & new_n70;
  assign new_n76 = new_n48 | x29 | (x4 & y0);
  assign new_n77 = new_n76 | (y1 & new_n62);
  assign new_n78 = new_n49 | x28 | (x2 & y0);
  assign new_n79 = new_n78 | (y1 & new_n63);
  assign new_n80 = new_n50 | x27 | (x0 & y0);
  assign new_n81 = new_n80 | (y1 & new_n64);
  assign new_n82 = new_n46 | x34 | (new_n47 & y0);
  assign new_n83 = new_n82 | (y1 & new_n54);
  assign new_n84 = new_n44 | x35 | (new_n45 & y0);
  assign new_n85 = new_n84 | (y1 & new_n55);
  assign new_n86 = new_n83 & new_n85;
  assign new_n87 = new_n86 & new_n81 & new_n79 & new_n75 & new_n77;
  assign y2 = ~new_n87;
  assign new_n89 = new_n50 | (y0 & x0) | (x18 & y1);
  assign new_n90 = new_n89 | (y2 & x27);
  assign new_n91 = new_n49 | (y0 & x2) | (x19 & y1);
  assign new_n92 = new_n91 | (y2 & x28);
  assign new_n93 = ~new_n90 | ~new_n92;
  assign new_n94 = new_n39 | (y0 & x8) | (x22 & y1);
  assign new_n95 = ~new_n94 & (~x31 | ~y2);
  assign new_n96 = ~new_n95;
  assign new_n97 = new_n41 | (y0 & x10) | (x23 & y1);
  assign new_n98 = new_n97 | (y2 & x32);
  assign new_n99 = new_n96 & new_n98;
  assign new_n100 = new_n40 | (y0 & x12) | (x24 & y1);
  assign new_n101 = ~new_n100 & (~x33 | ~y2);
  assign new_n102 = ~new_n101;
  assign new_n103 = new_n46 | (y0 & x14) | (x25 & y1);
  assign new_n104 = new_n103 | (y2 & x34);
  assign y4 = ~new_n104 | ~new_n99 | ~new_n102;
  assign new_n106 = new_n38 | (y0 & x6) | (x21 & y1);
  assign new_n107 = new_n106 | (y2 & x30);
  assign new_n108 = new_n48 | (y0 & x4) | (x20 & y1);
  assign new_n109 = ~new_n108 & (~x29 | ~y2);
  assign new_n110 = ~new_n109;
  assign new_n111 = ~new_n107 | ~new_n110;
  assign new_n112 = new_n44 | (y0 & x16) | (x26 & y1);
  assign new_n113 = new_n112 | (y2 & x35);
  assign y3 = new_n113 & (new_n93 | y4 | new_n111);
  assign new_n115 = ~new_n95 & ~new_n106 & (~x30 | ~y2);
  assign new_n116 = ~new_n98 | ~new_n102;
  assign new_n117 = ~new_n116;
  assign new_n118 = ~new_n104 | ~new_n102;
  assign y5 = new_n118 | (new_n109 & new_n99) | (new_n115 & new_n117);
  assign new_n120 = new_n104 & (new_n98 | new_n101);
  assign new_n121 = new_n120 & (new_n92 | new_n109 | new_n116) & (new_n116 | new_n107 | new_n95);
  assign y6 = ~new_n121;
endmodule



//+
Field[1].CurvesList = {1};
//+
Field[1].PointsList = {1};
//+
Field[1].CurvesList = {2};
//+
Field[1].PointsList = {1, 2};
//+
Field[2] = Constant;
//+
Delete Field [2];
//+
Field[1].PointsList = {2, 3};
//+
Delete Field [1];
//+
Field[1] = Distance;
//+
Field[1].CurvesList = {2};
//+
Field[1].PointsList = {2, 3};

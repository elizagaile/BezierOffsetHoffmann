using System;
using System.Collections.Generic;
using System.Drawing;
using System.Windows.Forms;
using MathNet.Numerics.LinearAlgebra; // matricām un vektoriem
// using MathNet.Numerics.RootFinding; // nepieciešams funkcijai "FindXfromt"

namespace OffsetHoffmann1
{
    public partial class OffsetForm : Form
    {
        Matrix<double> P = Matrix<double>.Build.Dense(4, 2); // matrica, kas satur galv enās līknes kontrolpunktus

        double maxError = 0.5; // vīles attāluma nepieciešamā precizitāte, te norādīta pikseļos
        double distance = -18; // ja negatīvs - kreisajā pusē no pirmā kotrolpunkta, pozitīvs - pa labi
        int movingPoint = -1; // nepieciešams līknes pārvietošanai ar peli
        List<double> splitValues = new List<double>(); // t vērtības, kurās nepieciešams sadalīt bezjē līkni

        // inicalizācija
        public OffsetForm()
        {
            InitializeComponent();
            
            // sākuma kontrolpunktu koordinātas:
            P[0, 0] = 83; 
            P[0, 1] = 303;

            P[1, 0] = 105;
            P[1, 1] = 162;

            P[2, 0] = 211;
            P[2, 1] = 186;

            P[3, 0] = 300;
            P[3, 1] = 309;

            // sākumā bezjē līkni sadala pa galapunktiem un ekstrēmiem
            splitValues.Add(0);
            splitValues.AddRange(BezierExtrema(P));
            splitValues.Add(1);
        }

        // attēlo visu grafiski, izsauc nepieciešamās funkcijas
        private void Canva_Paint(object sender, PaintEventArgs e)
        {
            e.Graphics.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias; // grafiski attēlo līnijas gludākas
            int r = 2; // punktu rādiuss

            // uzzīmē oriģinālo līkni, kontrolpunktus, poligonu
            for (int i = 0; i < MatrixToPointArray(P).Length; i++)
            {
                e.Graphics.DrawEllipse(Pens.Black, MatrixToPointArray(P)[i].X - r, MatrixToPointArray(P)[i].Y - r, 2 * r, 2 * r);
            }
            e.Graphics.DrawLines(Pens.LightGray, MatrixToPointArray(P));
            e.Graphics.DrawBeziers(Pens.Black, MatrixToPointArray(P));
            
            while (!IsSafe()) { } // atkārtoti izsauc funkciju, līdz līkne sadalīta tādos gabalos, ka sasniegta nepieciešamā precizitāte
            
            // sadales vērtības sakārto, un katram segmentam uzzīmē paralēli
            splitValues.Sort();
            for (int i = 0; i < splitValues.Count - 1; i++)
            {
                Matrix<double> offset = HoffmannOffset(SplitBezier(splitValues[i], splitValues[i + 1]));
                e.Graphics.DrawBeziers(Pens.Red, MatrixToPointArray(offset));
            }
            
            error.Text = "" + (splitValues.Count - 1); // lai redzētu, cik segmentos līkne sadalīta
        }

        // aprēķina paralēles kontrolpunktus, algoritms aprakstīts http://docs-hoffmann.de/bezier18122002.pdf
        private Matrix<double> HoffmannOffset(Matrix<double> P)
        {
            var M = Matrix<double>.Build;
            var V = Vector<double>.Build;

            var P0 = P.Row(0);
            var P1 = P.Row(1);
            var P2 = P.Row(2);
            var P3 = P.Row(3);

            var s0 = P1 - P0;
            var s3 = P3 - P2;

            // apskata, vai līkne nav deģenerēta - vai kādi punkti sakrīt
            if (P0.Equals(P1))
            {
                if (!P0.Equals(P2))
                {
                    s0 = P2 - P0;
                }
                else if (!P0.Equals(P3))
                {
                    s0 = P3 - P0;
                }
                else // visi punkti sakrīt
                {
                    MessageBox.Show("Degenerated Bezier curve does not have an offset!" + P);
                    return P;
                }
            }
            // apskata, vai līkne nav deģenerēta - vai kādi punkti sakrīt
            if (P3.Equals(P2))
            {
                if (!P3.Equals(P1))
                {
                    s3 = P3 - P1;
                }
                else if (!P3.Equals(P0))
                {
                    s3 = P3 - P0;
                }
                else // visi punkti sakrīt
                {
                    MessageBox.Show("Degenerated Bezier curve does not have an offset!");
                    return P;
                }
            }

            var a = 3 * s0 + 3 * s3 - 2 * (P3 - P0);
            var b = -6 * s0 - 3 * s3 + 3 * (P3 - P0);
            var c = 3 * s0;

            var n0 = V.Dense(2);
            n0[0] = -s0[1] / s0.L2Norm();
            n0[1] = s0[0] / s0.L2Norm();

            var n3 = V.Dense(2);
            n3[0] = -s3[1] / s3.L2Norm();
            n3[1] = s3[0] / s3.L2Norm();

            var Q0 = P0 + distance * n0;
            var Q1 = P1 + distance * n0;
            var Q2 = P2 + distance * n3;
            var Q3 = P3 + distance * n3;

            var A = 3 * s0 + 3 * s3 - 2 * (Q3 - Q0);
            var B = -6 * s0 - 3 * s3 + 3 * (Q3 - Q0);
            var C = 3 * s0;

            var Pc = 0.125 * a + 0.25 * b + 0.5 * c + P0;
            var Qc = 0.125 * A + 0.25 * B + 0.5 * C + Q0;

            var dP = 0.75 * a + b + c;
            var dQ = 0.75 * A + B + C;

            var nc = V.Dense(2);
            if (dP.L2Norm() != 0) // deģenerētā gadījumā var gadīties, ka līknes pieskare ir nulvektors
            {
                nc[0] = -dP[1] / dP.L2Norm();
                nc[1] = dP[0] / dP.L2Norm();
            }
            else
            {
                nc[0] = 0;
                nc[1] = 0;
            }

            var Rc = Pc + distance * nc;

            var tmp1 = s0.ToColumnMatrix().Transpose() * nc;
            var tmp2 = s3.ToColumnMatrix().Transpose() * nc;
            var tmp3 = 4 * (s0 - s3).ToRowMatrix() * nc;

            double[,] tmp = {{ s0[0], -s3[0], 8/3 * dP[0] },
                           { s0[1], -s3[1], 8/3 * dP[1] },
                           { tmp1[0], tmp2[0], tmp3[0] }};
            var matrixA = M.DenseOfArray(tmp);

            var y = V.Dense(3);
            y[0] = 8 / 3 * (Rc[0] - Qc[0]);
            y[1] = 8 / 3 * (Rc[1] - Qc[1]);
            tmp1 = 4 / 3 * dQ.ToColumnMatrix().Transpose() * nc;
            y[2] = tmp1[0];


            double dk0 = 0;
            double dk3 = 0;
            // vēl ģeometriski neatrisinātos gadījumos var gadīties, ka determinants ir nulle
            // tādā gadījumā pielieto naivu pārbīdi - dk0, dk3 = 0
            if (matrixA.Determinant() != 0)
            {
                var x = matrixA.Inverse() * y;
                dk0 = x[0];
                dk3 = x[1];
            }

            var R0 = Q0;
            var R1 = Q0 + (1 + dk0) * s0;
            var R2 = Q3 - (1 + dk3) * s3;
            var R3 = Q3;

            var R = M.DenseOfRowVectors(R0, R1, R2, R3);

            return R;
        }

        // pārbauda, vai paralēle ir pietikami precīza, ja nav, sadala smalkāk
        private bool IsSafe()
        {
            splitValues.Sort();
            double newSplitValue = -1;
            double checkpointCount = 20; // ātrdarbības nolūkos, paralēles precizitāti pārbauda diskrētā punktu skaitā

            for (int i = 0; i < splitValues.Count - 1; i++)
            {
                Matrix<double> segment = SplitBezier(splitValues[i], splitValues[i + 1]);
                Matrix<double> offset = HoffmannOffset(segment);
                
                for (double t = 0; t <= 1; t += 1 / checkpointCount)
                {
                    double length = GetLength(PointAtValue(t, segment), PointAtValue(t, offset));
                    if (Math.Abs(length - Math.Abs(distance)) > maxError)
                    {
                        newSplitValue = (splitValues[i] + splitValues[i + 1]) / 2; // ja paralēles segments nav pietiekmi precīzs, to sadala uz pusēm
                        break;
                    }
                }
            }

            if (newSplitValue < 0)
            {
                return true;
            }

            else
            {
                splitValues.Add(newSplitValue);
                return false;
            }
        }

        // atrod kontrolpunktu koordinātas bezjē līknes intervālam starp divam t vērtībām
        private Matrix<double> SplitBezier(double t1, double t2)
        {
            Matrix<double> res = Matrix<double>.Build.Dense(4, 2);

            res.SetRow(0, P.Row(0));
            res.SetRow(1, t2 * P.Row(1) - (t2 - 1) * P.Row(0));
            res.SetRow(2, t2 * t2 * P.Row(2) - 2 * t2 * (t2 - 1) * P.Row(1) + (t2 - 1) * (t2 - 1) * P.Row(0));
            res.SetRow(3, Math.Pow(t2, 3) * P.Row(3) - 3 * t2 * t2 * (t2 - 1) * P.Row(2) +
                          3 * t2 * (t2 - 1) * (t2 - 1) * P.Row(1) - Math.Pow(t2 - 1, 3) * P.Row(0));

            t1 = t1 / t2;

            res.SetRow(0, Math.Pow(t1, 3) * res.Row(3) - 3 * t1 * t1 * (t1 - 1) * res.Row(2) +
              3 * t1 * (t1 - 1) * (t1 - 1) * res.Row(1) - Math.Pow(t1 - 1, 3) * res.Row(0));
            res.SetRow(1, t1 * t1 * res.Row(3) - 2 * t1 * (t1 - 1) * res.Row(2) + (t1 - 1) * (t1 - 1) * res.Row(1));
            res.SetRow(2, t1 * res.Row(3) - (t1 - 1) * res.Row(2));
            res.SetRow(3, res.Row(3));

            return res;
        }

        // atrod bezjē līknes ekstrēmus
        private List<double> BezierExtrema(Matrix<double> M)
        {
            List<double> tValues = new List<double>();

            PointF[] derivitive = BezierDerivative(M);

            PointF a = new PointF
            {
                X = derivitive[0].X - 2 * derivitive[1].X + derivitive[2].X,
                Y = derivitive[0].Y - 2 * derivitive[1].Y + derivitive[2].Y
            };
            PointF b = new PointF
            {
                X = 2 * (derivitive[1].X - derivitive[0].X),
                Y = 2 * (derivitive[1].Y - derivitive[0].Y)
            };
            PointF c = new PointF
            {
                X = derivitive[0].X,
                Y = derivitive[0].Y
            };

            double t1 = (-b.X + Math.Sqrt(b.X * b.X - 4 * a.X * c.X)) / (2 * a.X);
            if (t1 > 0 && t1 < 1)
            {
                tValues.Add(t1);
            }

            double t2 = (-b.Y + Math.Sqrt(b.Y * b.Y - 4 * a.Y * c.Y)) / (2 * a.Y);
            if (t2 > 0 && t2 < 1)
            {
                tValues.Add(t2);
            }

            double t3 = (-b.X - Math.Sqrt(b.X * b.X - 4 * a.X * c.X)) / (2 * a.X);
            if (t3 > 0 && t3 < 1)
            {
                tValues.Add(t3);
            }

            double t4 = (-b.Y - Math.Sqrt(b.Y * b.Y - 4 * a.Y * c.Y)) / (2 * a.Y);
            if (t4 > 0 && t4 < 1)
            {
                tValues.Add(t4);
            }

            tValues.Sort();

            //pārbauda, vai kādi ekstrēmi nesakrīt; ja sakrīt, liekos izņem ārā
            for (int i = 0; i < tValues.Count - 1; i++)
            {
                if (tValues[i] == tValues[i + 1])
                {
                    tValues.RemoveAt(i + 1);
                    i--;
                }
            }

            return tValues;
        }

        // atrod bezjē līknes atvasinājuma koeficientus
        private PointF[] BezierDerivative(Matrix<double> M)
        {
            PointF[] derivitive = new PointF[3];

            for (int i = 0; i < M.RowCount - 1; i++)
            {
                PointF tmp = new PointF
                {
                    X = 3 * (float)(M[i + 1, 0] - M[i, 0]),
                    Y = 3 * (float)(M[i + 1, 1] - M[i, 1])
                };
                derivitive[i] = tmp;
            }

            return derivitive;
        }
        
        // aprēķina attālumu starp diviem punktiem
        private float GetLength(PointF firstPoint, PointF secondPoint)
        {
            return (float)Math.Sqrt(Math.Pow(firstPoint.X - secondPoint.X, 2) + Math.Pow(firstPoint.Y - secondPoint.Y, 2));
        }

        // atrod bezjē līknes koordinātas dotai t vērtībai
        private PointF PointAtValue(double t, Matrix<double> M)
        {
            double resultX = M[0, 0] * Math.Pow(1 - t, 3)
                           + 3 * M[1, 0] * Math.Pow(1 - t, 2) * t
                           + 3 * M[2, 0] * (1 - t) * t * t
                           + M[3, 0] * Math.Pow(t, 3);

            double resultY = M[0, 1] * Math.Pow(1 - t, 3)
                           + 3 * M[1, 1] * Math.Pow(1 - t, 2) * t
                           + 3 * M[2, 1] * (1 - t) * t * t
                           + M[3, 1] * Math.Pow(t, 3);

            return new PointF((float)resultX, (float)resultY);
        }

        // no n x 2 matricas izveido masīvu ar punktiem
        private PointF[] MatrixToPointArray(Matrix<double> M)
        {
            PointF[] points = new PointF[M.RowCount];

            if (M.ColumnCount == 2)
            {
                for (int i = 0; i < M.RowCount; i++)
                {
                    points[i] = new PointF { X = (float)M[i, 0], Y = (float)M[i, 1] };
                }
            }

            return points;
        }
        
        // atrod, vai un kurš kontrolpunkts ir peles kursora tuvumā
        private int FindLocalPoint(PointF mouseLocation)
        {
            const int localRadius = 7; // apkārtnes rādiusa izmērs pikseļos

            for (int i = 0; i < MatrixToPointArray(P).Length; i++)
            {
                if (GetLength(mouseLocation, MatrixToPointArray(P)[i]) < localRadius)
                {
                    return i;
                }
            }
            return -1;
        }

        // atrod, vai pele tika nospiesta tuvumā kontrolpunktam
        private void Canva_MouseDown(object sender, MouseEventArgs e)
        {
            movingPoint = FindLocalPoint(e.Location);
        }

        // ja izvēlēts kontrolpunkts, veic tā pārbīdi ar peli
        private void Canva_MouseMove(object sender, MouseEventArgs e)
        {
            if (movingPoint != -1)
            {
                P[movingPoint, 0] = e.X;
                P[movingPoint, 1] = e.Y;
                Canva.Invalidate();
            }
        }

        // pārtrauc kontrolpunkta pārbīdi
        private void Canva_MouseUp(object sender, MouseEventArgs e)
        {
            movingPoint = -1;
        }

        // ļāuj mainīt paralēles attālumu, izmantojot klavietūras uz augšu/leju bultiņas
        private void OffsetForm_KeyDown(object sender, KeyEventArgs e)
        {
            if (e.KeyCode == Keys.Up)
            {
                distance++;
                Canva.Invalidate();
            }

            else if (e.KeyCode == Keys.Down)
            {
                distance--;
                Canva.Invalidate();
            }
        }

        /*
        private Matrix<double> Adjoint(Matrix<double> P)
        {
            var M = Matrix<double>.Build;
            var C = M.Dense(3, 3);

            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    C[i, j] = P[(i + 2) % 3, (j + 2) % 3] * P[(i + 1) % 3, (j + 1) % 3] -
                    P[(i + 2) % 3, (j + 1) % 3] * P[(i + 1) % 3, (j + 2) % 3];
                }
            }

            return C.Transpose();
        }
        */
        
        /*
        private double MaxError(Matrix<double> offset)
        {
            double maxError = 0;
            for (double t = 0; t <= 1; t += 0.01)
            {
                List<double> tOffset = SegmentBezierIntersection(PointAtValue(t, P), NormalVector(t, P), offset);

                double correctError = -double.MaxValue;
                for (int i = 0; i < tOffset.Count; i++)
                {
                    correctError = Math.Min(Math.Abs(correctError), Math.Abs(GetLength(PointAtValue(t, P), PointAtValue(tOffset[i], offset)) - distance));
                }

                maxError = Math.Max(maxError, correctError);
            }
            return maxError;
        }
        */

        /*
        private List<double> SegmentBezierIntersection(PointF o, Vector<double> v, Matrix<double> M)
        {
            double alpha = 0;
            Matrix<double> Pnew = M.Clone();

            if (v[0] != 0)
            {
                alpha = Math.PI / 2 - Math.Atan(v[1] / v[0]);
            }
            
            for (int i = 0; i < M.RowCount; i++)
            {
                PointF tmp = new PointF()
                {
                    X = (float)M[i, 0],
                    Y = (float)M[i, 1]
                };
                Pnew[i, 0] = PointRotation(alpha, tmp, o).X;
                Pnew[i, 1] = PointRotation(alpha, tmp, o).Y;
            }

            Tuple<double, double, double> rootsTuple = FindtFromX(o.X, Pnew);
            List<double> roots = new List<double>
            {
                rootsTuple.Item1,
                rootsTuple.Item2,
                rootsTuple.Item3
            };
            for (int i = roots.Count - 1; i >= 0; i--)
            {
                if (Double.IsNaN(roots[i]) == true || roots[i] < 0 || roots [i] > 1)
                {
                    roots.RemoveAt(i);
                }
            }

            return roots;
        }
        */

        /*
        private Vector<double> NormalVector(double t, Matrix<double> M)
        {
            PointF old = PointAtValue(t, M);
            PointF[] derivCoeff = BezierDerivative(M);
            PointF derivitive = new PointF
            {
                X = (float)(derivCoeff[0].X * (1 - t) * (1 - t) + derivCoeff[1].X * 2 * (1 - t) * t + derivCoeff[2].X * t * t),
                Y = (float)(derivCoeff[0].Y * (1 - t) * (1 - t) + derivCoeff[1].Y * 2 * (1 - t) * t + derivCoeff[2].Y * t * t)
            };

            double length = Math.Sqrt(derivitive.X * derivitive.X + derivitive.Y * derivitive.Y);
            
            Vector<double> n = Vector<double>.Build.Dense(2);
            n[0] = (-derivitive.Y / length);
            n[1] = (derivitive.X / length);
            
            return n;
        }
        */
        
        /*
        private Tuple<double, double, double> FindtFromX(double x, Matrix<double> M)
        {
            double a = M[0, 0];
            double b = M[1, 0];
            double c = M[2, 0];
            double d = M[3, 0];

            double k3 = -a + 3 * b - 3 * c + d;
            double k2 = (3 * a - 6 * b + 3 * c) / k3;
            double k1 = (-3 * a + 3 * b) / k3;
            double k0 = (a - x) / k3;

            Tuple<double, double, double> roots = Cubic.RealRoots(k0, k1, k2);

            return roots;
        }
        */

        /*
        private PointF PointRotation(double alpha, PointF point, PointF pivot)
        {
            point.X -= pivot.X;
            point.Y -= pivot.Y;
            
            double xnew = point.X * Math.Cos(alpha) - point.Y * Math.Sin(alpha);
            double ynew = point.X * Math.Sin(alpha) + point.Y * Math.Cos(alpha);
            
            point.X = (float)xnew + pivot.X;
            point.Y = (float)ynew + pivot.Y;

            return point;
        }
        */
    }
}

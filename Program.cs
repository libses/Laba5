using ScottPlot;

namespace ConsoleApp40
{
    internal class Program
    {
        static void Main(string[] args)
        {
            var toCalc = new Func<double, double, double>((x, y) => y * y + x);
            var start = -2;
            CalculateAllEuler(toCalc, start);
            CalculateAllAdams(toCalc, start);
        }

        private static void CalculateAdams(Func<double, double, double> f, double start, int steps, System.Drawing.Color color, Plot plt)
        {
            var adamsSolver = new AdamsThreeSteps();
            var resAdams = adamsSolver.Solve(f, start, steps);
            (var dataAdams, var dataAdamsY) = resAdams.Convert();
            plt.AddScatter(dataAdams, dataAdamsY, color);
        }

        private static void CalculateAllAdams(Func<double, double, double> f, double start)
        {
            var plt = new Plot(1920, 1080);
            CalculateAdams(f, start, 10, System.Drawing.Color.Red, plt);
            CalculateAdams(f, start, 20, System.Drawing.Color.Green, plt);
            CalculateAdams(f, start, 30, System.Drawing.Color.Blue, plt);
            plt.SaveFig("adams.png");
        }

        private static void CalculateEuler(Func<double, double, double> f, double start, int steps, System.Drawing.Color color, Plot plt)
        {
            var eulerSolver = new EulerYavniy();
            var resEuler = eulerSolver.Solve(f, start, steps);
            (var dataEulerX, var dataEulerY) = resEuler.Convert();
            plt.AddScatter(dataEulerX, dataEulerY, color);
        }

        private static void CalculateAllEuler(Func<double, double, double> f, double start)
        {
            var plt = new Plot(1920, 1080);
            CalculateEuler(f, start, 10, System.Drawing.Color.Red, plt);
            CalculateEuler(f, start, 20, System.Drawing.Color.Green, plt);
            CalculateEuler(f, start, 30, System.Drawing.Color.Blue, plt);
            plt.SaveFig("euler.png");
        }
    }

    public static class Helper
    {
        public static (double[], double[]) Convert(this Point[] array)
        {
            return (array.Select(x => x.X).ToArray(), array.Select(x => x.Y).ToArray());
        }
    }

    public class Point
    {
        public double X;
        public double Y;
    }

    public interface IKoshiTaskSolver
    {
        Point[] Solve(Func<double, double, double> derivative, double valueInZero, int steps);
    }

    public static class RungeKuttaSpeedUp
    {
        public static double GetYnp1(double xn, double yn, double h, Func<double, double, double> f)
        {
            return yn + (h / 6) *
                (GetK1(xn, yn, h, f) +
                2 * GetK2(xn, yn, h, f) +
                2 * GetK3(xn, yn, h, f) +
                GetK4(xn, yn, h, f));
        }

        private static double GetK1(double xn, double yn, double h, Func<double, double, double> f)
        {
            return f(xn, yn);
        }

        private static double GetK2(double xn, double yn, double h, Func<double, double, double> f)
        {
            return f(xn + h/2, yn + GetK1(xn, yn, h, f) * h / 2);
        }

        private static double GetK3(double xn, double yn, double h, Func<double, double, double> f)
        {
            return f(xn + h / 2, yn + GetK2(xn, yn, h, f) * h / 2);
        }

        private static double GetK4(double xn, double yn, double h, Func<double, double, double> f)
        {
            return f(xn + h, yn + h * GetK3(xn, yn, h, f));
        }
    }

    public class AdamsThreeSteps : IKoshiTaskSolver
    {
        public Point[] Solve(Func<double, double, double> derivative, double valueInZero, int steps)
        {
            var h = 1d / steps;
            var points = new Point[steps + 1];
            points[0] = new Point() { X = 0, Y = valueInZero };
            var xnp1 = h;
            var xnp2 = 2 * h;
            var ynp1 = RungeKuttaSpeedUp.GetYnp1(0, valueInZero, h, derivative);
            var ynp2 = RungeKuttaSpeedUp.GetYnp1(xnp1, ynp1, h, derivative);
            points[1] = new Point() { X = xnp1, Y = ynp1 };
            points[2] = new Point() { X = xnp2, Y = ynp2 };
            for (int i = 3; i <= steps; i++)
            {
                var yi = GetYnp3(derivative, points[i - 1].Y, points[i - 1].X, points[i - 2].X, points[i - 2].Y, points[i - 3].X, points[i - 3].Y, h);
                points[i] = new Point { X = points[i - 1].X + h, Y = yi };
            }

            return points;
        }

        private double GetYnp3(Func<double, double, double> derivative, double ynp2, double xnp2, double xnp1, double ynp1, double xn, double yn, double h)
        {
            return ynp2 + h * ((23d / 12d) * derivative(xnp2, ynp2) - (16d / 12d) * derivative(xnp1, ynp1) + (5d / 12d) * derivative(xn, yn));
        } 
    }

    public class EulerYavniy : IKoshiTaskSolver
    {
        public Point[] Solve(Func<double, double, double> derivative, double valueInZero, int steps)
        {
            var points = new Point[steps + 1];
            var step = 1d / steps;
            var yi_1 = valueInZero;
            var xi_1 = 0d;
            var xi = 0d;
            var yi = 0d;

            for (int i = 0; i < steps; i++)
            {
                points[i] = new Point { X = xi_1, Y = yi_1 };
                xi = xi_1 + step;
                yi = GetYi(yi_1, xi, xi_1, derivative);
                yi_1 = yi;
                xi_1 = xi;
            }

            points[steps] = new Point { X = xi, Y = yi };
            return points;
        }

        private double GetYi(double yi_1, double xi, double xi_1, Func<double, double, double> derivative)
        {
            return yi_1 + (xi - xi_1) * derivative(xi_1, yi_1);
        }
    }
}
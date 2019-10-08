using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeneticAlgorithm
{
    class Start
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Start....");
            Ga genAlg = new Ga(30, 48, 1000, 0.8f, 0.9f);
            genAlg.Init("C:\\Users\\chenia\\Desktop\\data.txt");
            genAlg.Solve();
            Console.ReadKey();
        }
    }
}

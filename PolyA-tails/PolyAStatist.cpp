#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <cmath>
#include <algorithm>
//structure, representing repeat in genome
struct Repeat
{
  Repeat(std::string _family, int _start, int _end, std::string classOfRepeat, std::string chromName, char _strand)
  : family(_family), start(_start), end(_end), name(classOfRepeat), targetChromosome(chromName), strand(_strand), isTSD(false)
  { }
  Repeat()
    : family(""), start(0), end(0)
  { }

  std::string targetChromosome;
  std::string family;  //L1HS
  std::string name;    //LINE1
  int start;
  int end;
  char strand;
  std::string leftTSD;
  int startPosLeftTSD;
  std::string rightTSD;
  int startPosRightTSD;
  int lenghtOfTSD;
  int editDistance;
  bool isTSD;
  //print Repeat in format presented in https://github.com/oyasnev/tsddoc 
  void Print(std::ofstream& out)
  {
    out << targetChromosome << ";" << name << ";" << family << ";" << start << ";" << end << ";" << strand << ";";
    if (isTSD)
    {
      out << 1 << ";" << leftTSD << ";" << startPosLeftTSD << ";" << rightTSD << ";" << startPosRightTSD << ";" << lenghtOfTSD << ";" << editDistance << std::endl;
    }
    else
    {
      out << 0 << std::endl;
    }
  }
};
//read DNA from RepeatMasker .fa file
void readGenome(std::string& genome, std::ifstream& in)
{
  std::string temp;
  while (in >> temp)
  {
    genome.append(temp);
  }
  std::transform(genome.begin(), genome.end(), genome.begin(), ::toupper);
}

std::unordered_map<long long, Repeat> mapOfAnnotation;
//Read annotation of repeats from RepeatMasker fa.out file.
void fillAnnotation(std::string& filename)
{
  Repeat last;
  int numberOfpolyA = 0;
  int totalNumber = 0;
  std::ifstream in(filename);
  std::string temp = "";
  std::string temp2 = "";
  std::getline(in, temp);
  std::getline(in, temp);
  int startPos = 0;
  int endPos = 0;
  std::string familyOfRepeat = "";
  while (in >> temp)
  {
    for (int i = 0; i < 3; ++i)
      in >> temp;
    std::string chromName = "";
    in >> chromName;
    in >> startPos;
    in >> endPos;
    in >> temp;
    char strand;
    in >> strand;
    in >> familyOfRepeat;
    std::string classOfRepeat = "";
    in >> classOfRepeat;
    if (classOfRepeat.substr(0, 4) == "LINE")
    {
      Repeat newRep = Repeat(familyOfRepeat, startPos, endPos, classOfRepeat, chromName, strand);
      if (newRep.family == last.family && newRep.start - last.end < 10)
      {
        mapOfAnnotation[last.start].end = newRep.end;
        mapOfAnnotation[last.start].strand = newRep.strand;
      }
      else
      {
        mapOfAnnotation[startPos] = newRep;
      }
      last = newRep;
    }
    std::getline(in, temp);
  }
}


bool polyASearching(Repeat rep, std::string& genome)
{
  if (rep.strand == '+')
  {
    int startOfSearch = rep.end;
    while (startOfSearch >= 0 && (genome[startOfSearch] == 'A' || genome[startOfSearch] == 'T'))
    {
      startOfSearch--;
    }
    startOfSearch--;
    while (startOfSearch >= 0 && (genome[startOfSearch] == 'A' || genome[startOfSearch] == 'T'))
      startOfSearch--;
    int endOfSearch = rep.end + 1;
    while (endOfSearch < genome.length() && (genome[endOfSearch] == 'A' || genome[endOfSearch] == 'T'))
      endOfSearch++;
    endOfSearch++;
    while (endOfSearch < genome.length() && (genome[endOfSearch] == 'A' || genome[endOfSearch] == 'T'))
      endOfSearch++;
    if (std::abs(endOfSearch - startOfSearch) > 16)
    {
      return true;
    }
  }
  else
  {
    int startOfSearch = rep.start;
    while (startOfSearch >= 0 && (genome[startOfSearch] == 'A' || genome[startOfSearch] == 'T'))
      startOfSearch--;
    startOfSearch--;
    while (startOfSearch >= 0 && (genome[startOfSearch] == 'A' || genome[startOfSearch] == 'T'))
      startOfSearch--;
    int endOfSearch = rep.start + 1;
    while (endOfSearch < genome.length() && (genome[endOfSearch] == 'A' || genome[endOfSearch] == 'T'))
      endOfSearch++;
    endOfSearch++;
    while (endOfSearch < genome.length() && (genome[endOfSearch] == 'A' || genome[endOfSearch] == 'T'))
      endOfSearch++;
    if (std::abs(endOfSearch - startOfSearch) > 16)
    {
      return true;
    }
  }
  return false;

}

//Checking repeat for long enough poly-A tail
bool havePolyATail(Repeat rep, std::string& genome)
{
  //Sometimes we couldn't find poly-A tail and it could move a bit further form the end of repeat
  if (polyASearching(rep, genome))
    return true;

  if (rep.strand == '+')
  {
    for (int i = 0; i < 8; ++i)
    {
      rep.end += 20;
      if (polyASearching(rep, genome))
        return true;
    }
  }
  else
  {
    for (int i = 0; i < 8; ++i)
    {
      rep.start -= 20;
      if (polyASearching(rep, genome))
        return true;
    }
  }

  return false;
}

std::unordered_map<std::string, int> total;
std::unordered_map<std::string, int> polyAMap;

//Total statisctic
void findPolyA(std::string& genome)
{
  std::ofstream myOut("polACurr.txt");
  for (auto it : mapOfAnnotation)
  {
    //cut too short repeats, because they likely doesn't have poly-A tail
    if (it.second.end - it.second.start > 1000)
    {
      it.second.Print(myOut);
      if (havePolyATail(it.second, genome))
      {
        myOut << true << std::endl;
        polyAMap[it.second.family]++;
      }
      else
      {
        myOut << false << std::endl;
      }
      total[it.second.family]++;
    }
  }
}

void printPolyAStat(std::ofstream& out2, std::string& nameOfChrom, int num)
{
  if (num > 1)
  {
    out2 << "In " << nameOfChrom << " and " << num - 1 << " other found:" << std::endl;
  }
  else
  {
    out2 << "In " << nameOfChrom << " found:" << std::endl;
  }
  out2 << "Family Of Repeat;Number of repeats;Number of repeats with poly-A tails;Percentage of repeat with poly-A tails" << std::endl;
  for (auto it : total)
  {
    out2 << it.first << ";" << it.second << ";" << polyAMap[it.first] << ";" << (double)polyAMap[it.first] / (double)it.second << std::endl;
  }
}

int main(int argc, char* argv[])
{
  if (argc % 2 != 1)
  {
    std::cout << "Wrong args";
    return 0;
  }
  std::string nameOfChrom = "";
  for (int i = 0; i < argc - 1; i = i + 2)
  {
    std::string filenameOfAnot = argv[1];
    std::string filenameOfGenome = argv[2];
    std::ifstream inGenome(filenameOfGenome);
    std::string genome = "";
    inGenome >> nameOfChrom;
    readGenome(genome, inGenome);
    fillAnnotation(filenameOfAnot);
    findPolyA(genome);
    genome.clear();
    mapOfAnnotation.clear();
  }

  std::ofstream out2("polyA.out");
  printPolyAStat(out2, nameOfChrom, argc/2);
  return 0;
}


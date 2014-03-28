// TSDFinder.cpp: определяет точку входа для консольного приложения.
//

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iomanip>
//The idea is to use this program on human chromosomes
//with masked repeats by RepearMasker
//Assuming that TSD are masked, but not all TSD are known,
//we can try to find new TSD at the ends of repeat regions.
//Let's try!



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
  void Print(std::ofstream& out)
  {

    out << targetChromosome << ";" <<name << ";" << family << ";" << start << ";" << end << ";" << strand << ";";
    if(isTSD)
    {
      out << 1 << ";" << leftTSD << ";" << startPosLeftTSD << ";" << rightTSD << ";" << startPosRightTSD << ";" << lenghtOfTSD << ";" << editDistance << std::endl;
    }
    else
    {
      out << 0 << std::endl;
    }
  }
};



const int MAXPENALTY = 3;
const int INFTY = 1000;
std::unordered_map<long long, Repeat> mapOfAnnotation;
std::unordered_map<std::string, std::vector<std::string>> typesAndTSD;
std::unordered_map<std::string, int> expectedLengthOfRepeat;

void readLengthes(std::string filename)
{
  std::ifstream in(filename);
  std::string temp = "";
  int temp2 = 0;
  std::string temp3 = "";

  while(in >> temp)
  {
    in >>temp2 >> temp3;
    expectedLengthOfRepeat[temp3] = std::max(temp2, expectedLengthOfRepeat[temp3]);
  }
}

void fillAnnotation(std::string& filename)
{
  std::ifstream in(filename);
  std::string temp = "";
  std::string temp2 = "";
  std::getline(in, temp);
  std::getline(in, temp);
  int startPos = 0;
  int endPos = 0;
  std::string familyOfRepeat = "";
  while(in >> temp)
  {
    for(int i = 0; i < 3; ++i)
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
    if(classOfRepeat.substr(0, 4) == "LINE")
      mapOfAnnotation[startPos] = Repeat(familyOfRepeat, startPos, endPos, classOfRepeat, chromName, strand);
    std::getline(in, temp);
  }
}

void readGenome(std::string& genome, std::ifstream& in)
{
  std::string temp;
  while(in >> temp)
  {
    genome.append(temp);
  }
}


void alignOverlapTSD(std::string left, std::string right, std::ofstream& out, size_t startOfRepeat, Repeat myRepeat)
{
  if(std::min(left.length(), right.length()) < 6)
    return;
  std::vector<std::vector<int>> matrix(left.length() + 1, std::vector<int>(right.length() + 1));
  std::vector<std::vector<int>> path(left.length() + 1, std::vector<int>(right.length() + 1));
  int maxScore = -INFTY;
  int posX = 0;
  int posY = 0;
  for(int i = 1; i <= left.length(); ++i)
  {
    for(int j = 1; j <= right.length(); ++j)
    {
      if(left[i - 1] == right[j - 1])
      {
        matrix[i][j] = std::max(matrix[i - 1][j - 1] + 1, std::max(matrix[i - 1][j] - 3, std::max(matrix[i][j - 1] - 3, 0)));
        if(matrix[i][j] > maxScore)
        {
          maxScore = matrix[i][j];
          posX = i;
          posY = j;
        }
      }
      else
      {
        matrix[i][j] = std::max(matrix[i - 1][j - 1] - 3, std::max(matrix[i - 1][j] - 3, std::max(matrix[i][j - 1] - 3, 0)));
      }
    }
  }


  //Backtrack

  std::string answer1 = "";
  std::string answer2 = "";
  int xPos = posX;
  int yPos = posY;
  int numberOfLetterOfA = 0;
  int numberOfLetterOfB = 0;


  while(matrix[xPos][yPos] != 0)
  {
    if(left[xPos - 1] == right[yPos - 1] && matrix[xPos][yPos] == matrix[xPos - 1][yPos - 1] + 1)
    {
      answer1.push_back(left[xPos - 1]);
      answer2.push_back(right[yPos - 1]);
      xPos--;
      yPos--;
      numberOfLetterOfA++;
      numberOfLetterOfB++;
      continue;
    }

    if(left[xPos - 1] != right[yPos - 1] && matrix[xPos][yPos] == matrix[xPos - 1][yPos - 1] - 3)
    {
      answer1.push_back(left[xPos - 1]);
      answer2.push_back(right[yPos - 1]);
      xPos--;
      yPos--;
      numberOfLetterOfA++;
      numberOfLetterOfB++;
      continue;
    }

    if(matrix[xPos][yPos] == matrix[xPos - 1][yPos] - 3)
    {
      answer1.push_back(left[xPos - 1]);
      answer2.push_back('-');
      xPos--;
      numberOfLetterOfA++;
      continue;
    }

    if(matrix[xPos][yPos] == matrix[xPos][yPos - 1] - 3)
    {
      answer1.push_back('-');
      answer2.push_back(right[yPos - 1]);
      yPos--;
      numberOfLetterOfB++;
      continue;
    }
  }
  std::reverse(answer1.begin(), answer1.end());
  std::reverse(answer2.begin(), answer2.end());

  int distance = 0;
  for(int i = 0; i < answer1.length(); ++i)
  {
    if(answer1[i] != answer2[i])
    {
      distance++;
    }
  }

  if((!(distance >= 3)) && std::min(numberOfLetterOfA, numberOfLetterOfB) >= 5)
  {
    if(!(std::min(numberOfLetterOfA, numberOfLetterOfB) <= 6 && distance >= 2))
    {
      myRepeat.isTSD = true;
      myRepeat.leftTSD = answer1;
      myRepeat.rightTSD = answer2;
      myRepeat.editDistance = distance;
      myRepeat.lenghtOfTSD = answer1.length();
      myRepeat.startPosLeftTSD = myRepeat.start - left.length() + xPos;
      myRepeat.startPosRightTSD = myRepeat.end + yPos;
    }
  }

  myRepeat.Print(out);
  if(mapOfAnnotation.find(startOfRepeat + 1) != mapOfAnnotation.end())
  {
    typesAndTSD[mapOfAnnotation[startOfRepeat + 1].family].push_back(answer1 + " " + answer2);
  }

}


void findTSD(std::string filename)
{
  std::ifstream in(filename);
  std::ofstream out(filename + ".newout");

  std::string nameOfChrom;
  in >> nameOfChrom;
  std::string genome;
  readGenome(genome, in);

  for(auto it : mapOfAnnotation)
  {
    size_t startOfRegion = it.second.start - 1;
    size_t endOfRegion = it.second.end;
    size_t length = endOfRegion - startOfRegion;

    if((it.second.family == "L1HS" && length <= 6000) || it.second.name != "LINE/L1")
    {
      continue;
    }

    if(it.second.family != "L1HS" && ((expectedLengthOfRepeat[it.second.family] == 0 && length > 1000)  || (double)expectedLengthOfRepeat[it.second.family]/(double)length < 0.95))
    {
 //     std::cout << length << " "<< expectedLengthOfRepeat[it.second.family] << " " << it.second.family << std::endl;
      continue;
    }

    //As we know the usual TSD length is 6-20 bp. So we will find such a strings
    int leftBorderOfLeftTSD = std::max(0, (int)startOfRegion - 300);
    std::string leftTSD = genome.substr(leftBorderOfLeftTSD, startOfRegion - leftBorderOfLeftTSD);
    int rightBorderOfLeftTSD = std::min(genome.length() - 1, endOfRegion + 300);
    std::string rightTSD = genome.substr(endOfRegion, rightBorderOfLeftTSD - endOfRegion);

    std::transform(leftTSD.begin(), leftTSD.end(), leftTSD.begin(), ::toupper);
    std::transform(rightTSD.begin(), rightTSD.end(), rightTSD.begin(), ::toupper);    
    if(leftTSD != "")
    {
      if(it.second.family == "L1HS" )
        std::cout << "Hi! Trying to align myself - " << leftTSD << " and " << rightTSD << std::endl;

      alignOverlapTSD(leftTSD, rightTSD, out, startOfRegion, it.second);
    }
  }
  std::ofstream secondout(filename + ".second.out");
  for(auto it : typesAndTSD)
  {
    secondout << "--------------------" << it.first << "------------------------" << std::endl;
    std::unordered_map<std::string, int> tempMap;
    for(size_t i = 0; i < it.second.size(); ++i)
    {
      tempMap[it.second[i]]++;;
    }
    for(auto it2 : tempMap)
    {
      secondout << std::setw(20) << it2.first << "  " << it2.second << std::endl;
    }
  }
}


int main(int argc, char* argv[])
{
  if(argc != 4)
  {
    std::cout << "Wrong args";
    return 0;
  }
  std::string lengthFilename = argv[3];
  readLengthes(lengthFilename);
  std::string filename = argv[1];
  std::string secondFilename = argv[2];
  fillAnnotation(secondFilename);
  findTSD(filename);
	return 0;
}


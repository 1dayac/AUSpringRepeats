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
  Repeat(std::string _family, int _start, int _end)
    : family(_family), start(_start), end(_end)
  { }
  Repeat()
    : family(""), start(0), end(0)
  { }

  std::string family;
  int start;
  int end;
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
    expectedLengthOfRepeat[temp] = temp2;
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
    for(int i = 0; i < 4; ++i)
      in >> temp;
    in >> startPos;
    in >> endPos;
    for(int i = 0; i < 2; ++i)
      in >> temp;
    in >> familyOfRepeat;
    std::string classOfRepeat = "";
    in >> classOfRepeat;
    if(classOfRepeat.substr(0, 4) == "LINE" || classOfRepeat.substr(0, 4) == "SINE" || classOfRepeat.substr(0, 3) == "LTR")
      mapOfAnnotation[startPos] = Repeat(familyOfRepeat, startPos, endPos);
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


void alignOverlapTSD(std::string left, std::string right, std::ofstream& out, size_t startOfRepeat)
{
  if(std::min(left.length(), right.length()) < 6)
    return;
  std::vector<std::vector<int>> matrix(left.length() + 1, std::vector<int>(right.length() + 1));
  std::vector<std::vector<int>> path(left.length() + 1, std::vector<int>(right.length() + 1));
  for(int i = 1; i <= left.length(); ++i)
  {
    for(int j = 1; j <= right.length(); ++j)
    {
      if(left[i - 1] == right[j - 1])
      {
        matrix[i][j] = std::min(matrix[i - 1][j - 1], std::min(matrix[i - 1][j] + 1, matrix[i][j - 1] + 1));
      }
      else
      {
        matrix[i][j] = std::min(matrix[i - 1][j - 1] + 1, std::min(matrix[i - 1][j] + 1, matrix[i][j - 1] + 1));
      }
    }
  }

  int position = 0;
  int distance = INFTY;
  for(int i = 5; i < right.length(); ++i)
  {
    if(matrix[left.length()][i] <= distance)
    {
      position = i;
      distance = matrix[left.length()][i];
    }
  }

  if(distance > MAXPENALTY)
    return;
  //Backtrack
  std::string answer1 = "";
  std::string answer2 = "";
  int xPos = left.length();
  int yPos = position;
  int numberOfLetterOfA = 0;
  int numberOfLetterOfB = 0;

  while(xPos != 0 && yPos != 0)
  {
    if(left[xPos - 1] == right[yPos - 1] && matrix[xPos][yPos] == matrix[xPos - 1][yPos - 1])
    {
      answer1.push_back(left[xPos - 1]);
      answer2.push_back(right[yPos - 1]);
      xPos--;
      yPos--;
      numberOfLetterOfA++;
      numberOfLetterOfB++;
      continue;
    }

    if(left[xPos - 1] != right[yPos - 1] && matrix[xPos][yPos] == matrix[xPos - 1][yPos - 1] + 1)
    {
      answer1.push_back(left[xPos - 1]);
      answer2.push_back(right[yPos - 1]);
      xPos--;
      yPos--;
      numberOfLetterOfA++;
      numberOfLetterOfB++;
      continue;
    }

    if(matrix[xPos][yPos] == matrix[xPos - 1][yPos] + 1)
    {
      answer1.push_back(left[xPos - 1]);
 //     answer2.push_back('-');
      xPos--;
      numberOfLetterOfA++;
      continue;
    }

    if(matrix[xPos][yPos] == matrix[xPos][yPos - 1] + 1)
    {
 //     answer1.push_back('-');
      answer2.push_back(right[yPos - 1]);
      yPos--;
      numberOfLetterOfB++;
      continue;
    }
  }
  if(std::min(numberOfLetterOfA, numberOfLetterOfB) < 5)
    return;
  if(std::min(numberOfLetterOfA, numberOfLetterOfB) <= 6 && distance >= 2)
    return;
  std::reverse(answer1.begin(), answer1.end());
  std::reverse(answer2.begin(), answer2.end());
  if(mapOfAnnotation.find(startOfRepeat + 1) != mapOfAnnotation.end())
  {
    out << answer1 << " " << answer2 << " " << mapOfAnnotation[startOfRepeat + 1].family << std::endl;
    typesAndTSD[mapOfAnnotation[startOfRepeat + 1].family].push_back(answer1);
  }
}


void findTSD(std::string filename)
{
  std::ifstream in(filename);
  std::ofstream out(filename + ".out");

  std::string nameOfChrom;
  in >> nameOfChrom;
  std::string genome;
  readGenome(genome, in);

  for(auto it : mapOfAnnotation)
  {
    size_t startOfRegion = it.second.start - 1;
    size_t endOfRegion = it.second.end;
    size_t length = endOfRegion - startOfRegion;
    if((double)expectedLengthOfRepeat[it.second.family]/(double)length > 1.1 || (double)expectedLengthOfRepeat[it.second.family]/(double)length < 0.9)
      continue;
    //As we know the usual TSD length is 6-20 bp. So we will find such a strings
    int leftBorderOfLeftTSD = std::max(0, (int)std::max(genome.find_last_of("N", startOfRegion - 1) + 1, startOfRegion - 20));
    std::string leftTSD = genome.substr(leftBorderOfLeftTSD, startOfRegion - leftBorderOfLeftTSD);
    int rightBorderOfLeftTSD = std::min(genome.length() - 1, std::min(genome.find_first_of("N", endOfRegion), endOfRegion + 20));
    std::string rightTSD = genome.substr(endOfRegion, rightBorderOfLeftTSD - endOfRegion);
    if(leftTSD != "")
      alignOverlapTSD(leftTSD, rightTSD, out, startOfRegion);
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


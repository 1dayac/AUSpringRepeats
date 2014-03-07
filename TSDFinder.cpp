// TSDFinder.cpp: определяет точку входа для консольного приложения.
//

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
//The idea is to use this program on human chromosomes
//with masked repeats by RepearMasker
//Assuming that TSD are masked, but not all TSD are known,
//we can try to find new TSD at the ends of repeat regions.
//Let's try!

const int MAXPENALTY = 3;
const int INFTY = 1000;
std::unordered_map<long long, std::string> mapOfAnnotation;
std::unordered_map<std::string, std::vector<std::string>> typesAndTSD;


void fillAnnotation(std::string& filename)
{
  std::ifstream in(filename);
  std::string temp = "";
  std::string temp2 = "";
  std::getline(in, temp);
  std::getline(in, temp);
  int startPos = 0;
  std::string familyOfRepeat = "";
  while(in >> temp)
  {
    for(int i = 0; i < 4; ++i)
      in >> temp;
    in >> startPos;
    for(int i = 0; i < 3; ++i)
      in >> temp;
    in >> familyOfRepeat;
    std::string classOfRepeat = "";
    in >> classOfRepeat;
    if(classOfRepeat.substr(0, 4) == "LINE" || classOfRepeat.substr(0, 4) == "SINE" || classOfRepeat.substr(0, 3) == "LTR")
      mapOfAnnotation[startPos] = familyOfRepeat;
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
      answer2.push_back('-');
      xPos--;
      numberOfLetterOfA++;
      continue;
    }

    if(matrix[xPos][yPos] == matrix[xPos][yPos - 1] + 1)
    {
      answer1.push_back('-');
      answer2.push_back(right[yPos - 1]);
      yPos--;
      numberOfLetterOfB++;
      continue;
    }
  }
  if(std::min(numberOfLetterOfA, numberOfLetterOfB) < 5)
    return;
  std::reverse(answer1.begin(), answer1.end());
  std::reverse(answer2.begin(), answer2.end());
  if(mapOfAnnotation.find(startOfRepeat + 1) != mapOfAnnotation.end())
  {
    out << answer1 << " " << answer2 << " " << mapOfAnnotation[startOfRepeat + 1] << std::endl;
    typesAndTSD[mapOfAnnotation[startOfRepeat + 1]].push_back(answer1);
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
  for(size_t i = 0; i < genome.length(); ++i)
  {
    size_t startOfRegion = genome.find_first_of("N", i);
    size_t endOfRegion = genome.find_first_not_of("N", startOfRegion);
    if(endOfRegion == std::string::npos)
      break;
    i = endOfRegion;
    //As we know the usual TSD length is 6-20 bp. So we will find such a strings
    int leftBorderOfLeftTSD = std::max(0, (int)std::max(genome.find_last_of("N", startOfRegion - 1) + 1, startOfRegion - 20));
    std::string leftTSD = genome.substr(leftBorderOfLeftTSD, startOfRegion - leftBorderOfLeftTSD);
    int rightBorderOfLeftTSD = std::min(genome.length() - 1, std::min(genome.find_first_of("N", endOfRegion), endOfRegion + 20));
    std::string rightTSD = genome.substr(endOfRegion, rightBorderOfLeftTSD - endOfRegion);
    alignOverlapTSD(leftTSD, rightTSD, out, startOfRegion);
  }
  std::ofstream secondout(filename + ".second.out");
  for(auto it : typesAndTSD)
  {
    secondout << it.first << ": ";
    for(size_t i = 0; i < it.second.size(); ++i)
    {
      secondout << it.second[i] << " ";
    }
    secondout << std::endl;
  }
}


int main(int argc, char* argv[])
{
  if(argc != 3)
  {
    std::cout << "Wrong args";
    return 0;
  }

  std::string filename = argv[1];
  std::string secondFilename = argv[2];
  fillAnnotation(secondFilename);
  findTSD(filename);
	return 0;
}


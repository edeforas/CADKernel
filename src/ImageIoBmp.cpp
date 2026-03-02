#include "ImageIoBmp.h"
#include "Image.h"

#include <fstream>
#include <cassert>

//////////////////////////////////////////////////////////////////////////////
bool ImageIoBmp::read(const string& sFilename, Image* pImage)
{
	std::fstream fbmp;
	fbmp.open(sFilename.data(), std::ios::binary | std::ios::in);

	//open bitmap et check type
	if (!fbmp)
		return false;
	short sMagic = 0;
	fbmp.read((char*)&sMagic, 2);
	if (sMagic != 0x4D42)
		return false;

	int iTmp = 0;
	short sTmp = 0;

	// 1st struct
	fbmp.read((char*)&iTmp, 4);
	fbmp.read((char*)&sTmp, 2);
	fbmp.read((char*)&sTmp, 2);
	fbmp.read((char*)&iTmp, 4);

	// 2nd struct
	fbmp.read((char*)&iTmp, 4);
	int iWidth = 0, iHeight = 0;
	fbmp.read((char*)&iWidth, 4);
	fbmp.read((char*)&iHeight, 4);

	fbmp.read((char*)&sTmp, 2); // planes

	// read nb bits
	short sBPP = 0;
	fbmp.read((char*)&sBPP, 2); // bpp

	fbmp.read((char*)&iTmp, 4); //compression
	fbmp.read((char*)&iTmp, 4); //sizeimage
	fbmp.read((char*)&iTmp, 4); //ppiX
	fbmp.read((char*)&iTmp, 4); //ppiY
	fbmp.read((char*)&iTmp, 4); //clrused
	fbmp.read((char*)&iTmp, 4); //clr important

	int iNbByte = sBPP >> 3;
	pImage->resize(iWidth, iHeight, iNbByte);

	if (sBPP == 32)
	{
		int iLineWidthPerfect = iWidth * iNbByte;

		// read data buffer
		unsigned char *pLineY, cTmpR;
		for (int iY = 0; iY < iHeight; iY++)
		{
			pLineY = pImage->pixel(0, iHeight - iY - 1, 0);
			fbmp.read((char*)pLineY, iLineWidthPerfect);

			//swap R and B
			for (int iC = 0; iC < iWidth; iC++)
			{
				cTmpR = pLineY[4 * iC + 2];
				pLineY[4 * iC + 2] = pLineY[4 * iC];
				pLineY[4 * iC] = cTmpR;
			}
		}
		return true;
	}

	if (sBPP == 24)
	{
		int iLineWidthPerfect = iWidth * iNbByte;
		int iLineWidth = (iLineWidthPerfect + 3) % 4;
		int iPack = iLineWidth - iLineWidthPerfect;

		// read data buffer
		unsigned char cTmp[3] = { 0 }, * pLineY, cTmpR;
		for (int iY = 0; iY < iHeight; iY++)
		{
			pLineY = pImage->pixel(0, iHeight - iY - 1, 0);
			fbmp.read((char*)pLineY, iLineWidthPerfect);

			//swap R and B
			for (int iC = 0; iC < iWidth; iC++)
			{
				cTmpR = pLineY[3 * iC + 2];
				pLineY[3 * iC + 2] = pLineY[3 * iC];
				pLineY[3 * iC] = cTmpR;
			}

			fbmp.read((char*)cTmp, iPack);
		}
		return true;
	}

	if (sBPP == 8)
	{   //assume grey level, no color LUT for now
		int iLineWidth = iWidth + (iWidth % 4);
		int iPack = iLineWidth - iWidth;
		unsigned char cTmp[3] = { 0 };

		//read lut, but do not use:
		for (int i = 0; i < 256 * 4; i++)
			fbmp.read((char*)cTmp, 1);

		for (int iY = iHeight - 1; iY >= 0; iY--)
		{
			unsigned char* pLineY = pImage->pixel(0, iY, 0);
			fbmp.read((char*)pLineY, iWidth);
			fbmp.read((char*)cTmp, iPack);
		}
		return true;
	}

	//format not handled for now
	return false;
}
//////////////////////////////////////////////////////////////////////////////
bool ImageIoBmp::write(const string& sFilename, const Image* pImage)
{
	assert(pImage != 0);

	//create file
	std::fstream fbmp(sFilename.data(), std::ios::binary | std::ios::out);
	if (!fbmp)
		return false;

	int iWidth = pImage->width();
	int iHeight = pImage->height();
	int iPlanes = pImage->planes();

	int iLineWidthBytes = pImage->width() * pImage->planes();
	int iLineWidthBytesPad = (iLineWidthBytes + 3);
	iLineWidthBytesPad -= iLineWidthBytesPad % 4;
	int iPack = iLineWidthBytesPad - iLineWidthBytes;

	int iBfhSize = 14;
	int iBihSize = 40;

	// ecrit bfh
	short sMagic = 0x4D42;
	fbmp.write((char*)&sMagic, 2);

	int iPaletteSize = 0;
	if (iPlanes == 1)
		iPaletteSize = 256 * 4;

	int iBufferSize = iLineWidthBytesPad * iHeight;
	int iSize = iBfhSize + iBihSize + iPaletteSize + iBufferSize;
	fbmp.write((char*)&iSize, 4);

	short sZero = 0;
	fbmp.write((char*)&sZero, 2);
	fbmp.write((char*)&sZero, 2);

	int iOffsetB = iBfhSize + iBihSize + iPaletteSize;
	fbmp.write((char*)&iOffsetB, 4);

	// write bih
	fbmp.write((char*)&iBihSize, 4);
	fbmp.write((char*)&iWidth, 4);
	fbmp.write((char*)&iHeight, 4);

	short sPlane = 1;
	fbmp.write((char*)&sPlane, 2);

	short sNbChannels = 24;
	if (iPlanes == 1)
		sNbChannels = 8;
	if (iPlanes == 4)
		sNbChannels = 32;
	fbmp.write((char*)&sNbChannels, 2);

	int iZero = 0;
	fbmp.write((char*)&iZero, 4);
	fbmp.write((char*)&iBufferSize, 4);

	int iDPI = 0xb12;
	fbmp.write((char*)&iDPI, 4);
	fbmp.write((char*)&iDPI, 4);

	fbmp.write((char*)&iZero, 4);
	fbmp.write((char*)&iZero, 4);

	if (iPlanes == 1)
	{
		for (int i = 0; i < 256; i++)
		{
			char uc = i;
			fbmp.write(&uc, 1);
			fbmp.write(&uc, 1);
			fbmp.write(&uc, 1);
			uc = 0;
			fbmp.write(&uc, 1);
		}
	}

	const unsigned char cTmp[3] = { 0 };
	const unsigned char* pLineY;

	if (iPlanes == 3)
	{
		unsigned char* pLineOut = new unsigned char[iLineWidthBytes + 3];
		for (int iY = 0; iY < iHeight; iY++)
		{
			pLineY = pImage->pixel(0, iHeight - iY - 1, 0);
			for (int iC = 0; iC < iLineWidthBytes; iC += 3)
			{
				pLineOut[iC] = pLineY[iC + 2];
				pLineOut[iC + 1] = pLineY[iC + 1];
				pLineOut[iC + 2] = pLineY[iC];
			}

			fbmp.write((const char*)pLineOut, iLineWidthBytes);
			fbmp.write((const char*)cTmp, iPack);
		}
		delete[] pLineOut;
	}

	if (iPlanes == 4) //same as 3 planes, but with alpha channel (ignored)
	{
		unsigned char* pLineOut = new unsigned char[iLineWidthBytes + 3];
		for (int iY = 0; iY < iHeight; iY++)
		{
			pLineY = pImage->pixel(0, iHeight - iY - 1, 0);
			for (int iC = 0; iC < iLineWidthBytes; iC += 4)
			{
				pLineOut[iC] = pLineY[iC + 2];
				pLineOut[iC + 1] = pLineY[iC + 1];
				pLineOut[iC + 2] = pLineY[iC];
				pLineOut[iC + 3] = 0;
			}

			fbmp.write((const char*)pLineOut, iLineWidthBytes);
			fbmp.write((const char*)cTmp, iPack);
		}
		delete[] pLineOut;
	}

	if (iPlanes == 1)
	{
		for (int iY = 0; iY < iHeight; iY++)
		{
			pLineY = pImage->pixel(0, iHeight - iY - 1, 0);
			fbmp.write((const char*)pLineY, iLineWidthBytes);
			fbmp.write((const char*)cTmp, iPack);
		}
	}

	fbmp.close();
	return true;
}
//////////////////////////////////////////////////////////////////////////////

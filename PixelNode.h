#pragma once

class Cluster;

class PixelNode
{
public:
	PixelNode() { }

	PixelNode(float x, float y, float z, float g)
	{
		PixelNode::x = x;
		PixelNode::y = y;
		PixelNode::z = z;
		PixelNode::g = g;
	}

public:
	float x = -1.0f;
	float y = -1.0f;
	float z = -1.0f;
	float g = -1.0f;
	Cluster* parent = nullptr;
};
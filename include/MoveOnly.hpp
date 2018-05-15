#pragma once

class MoveOnly {
protected:
	bool valid = true;

	MoveOnly() {}

	MoveOnly(MoveOnly&&) = delete;
	MoveOnly(const MoveOnly&&) = delete;	
	MoveOnly(MoveOnly&) = delete;
	MoveOnly(const MoveOnly&) = delete;

	virtual ~MoveOnly() = 0;
};

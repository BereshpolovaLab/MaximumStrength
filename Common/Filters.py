#!/usr/bin/env python3
# 24, BereshpolovaLab, University of Connecticut


def MedianFilter(List):
    Value = 0
    if List:
        List.sort()
        Length = len(List)
        Value = List[int(Length / 2)] if Length % 2\
            else (List[int(Length / 2) - 1] + List[int(Length / 2)]) / 2
    return Value

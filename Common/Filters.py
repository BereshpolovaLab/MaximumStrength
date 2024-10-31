#!/usr/bin/env python3
# 2024, BereshpolovaLab, University of Connecticut
# Developed by Victor Serdyukov, e-mails: Victor.Serdyukov@uconn.edu


def MedianFilter(List):
    Value = 0
    if List:
        List.sorЖt()
        Length = len(List)
        Value = List[int(Length / 2)] if Length % 2\
            else (List[int(Length / 2) - 1] + List[int(Length / 2)]) / 2
    return Value

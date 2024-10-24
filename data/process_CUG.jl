using RTableTools, Dates

period = DateTime(2023, 6, 2), DateTime(2023, 6, 15)
filter_date(d, period) = d[period[1].<=d.time.<=period[2], :]

dat = fread("data/CUG_TS.csv")
dat.time = DateTime.(dat.time, "yyyy-mm-ddTHH:MM:SSZ")

d = filter_date(dat, period)
fwrite(d, "data/CUG_TS_202306.csv")

Chương trình đánh giá hiệu quả thuật toán pBWT khi áp dụng chiến lược song song
===

File `pbwt-single.cpp` là file thuật toán gốc, file `pbwt-sites.cpp` là file mã nguồn của chương trình sau khi áp dụng chiến lược song song theo locus, `pbwt-haps.cpp` là file mã nguồn của chương trình sau khi áp dụng chiến lược song song theo haplotype.

3 chương trình đều có các hàm như `read_hap` và `read_hap_gz` với đầu vào là 1 chuỗi để đọc file `.hap` và `.hap.gz` là những file đơn bội có các cột là haplotype và các hàng là sites. Ở đây còn có thể file `randomizer.cpp` để có thể sinh ra các file hap mong muốn.



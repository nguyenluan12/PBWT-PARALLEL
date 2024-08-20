Chương trình đánh giá hiệu quả thuật toán pBWT khi áp dụng chiến lược song song
===

File `pbwt_1.cpp` là file thuật toán gốc còn file `pbwt.cpp` là file mã nguồn của chương trình sau khi áp dụng chiến lược.

2 chương trình đều có các hàm như `read_hap` và `read_hap_gz` với đầu vào là 1 chuỗi để đọc file `.hap` và `.hap.gz` là những file đơn bội có các cột là haplotype và các hàng là sites. Ở đây còn có thể file `randomizer.cpp` để có thể sinh ra các file hap mong muốn.

Trong hàm main của 2 file `pbwt_1.cpp` và `pbwt.cpp` đã được làm sẵn để sửa tham số và đường dẫn là có thể sử dụng để đo đếm.
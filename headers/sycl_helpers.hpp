#pragma once
#include <memory>


namespace mdutil {

	template <typename Ty>
	inline decltype(auto) as_3d(const std::unique_ptr<Ty[]>& arr) {
		return reinterpret_cast<Ty(*)[3]>(arr.get());
	}
}
